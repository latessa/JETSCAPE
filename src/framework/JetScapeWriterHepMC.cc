/*******************************************************************************
 * Copyright (c) The JETSCAPE Collaboration, 2018
 *
 * Modular, task-based framework for simulating all aspects of heavy-ion
 *collisions
 *
 * For the list of contributors see AUTHORS.
 *
 * Report issues at https://github.com/JETSCAPE/JETSCAPE/issues
 *
 * or via email to bugs.jetscape@gmail.com
 *
 * Distributed under the GNU General Public License 3.0 (GPLv3 or later).
 * See COPYING for details.
 ******************************************************************************/

/*******************************************************************************
 * kk Sep 22, 2020:  Significant updates to status codes
 * to conform with Appendix A in https://arxiv.org/pdf/1912.08005.pdf
 * We will follow Pythia's example of preserving as much of the internal status
 code as possible
 * Specifically:
 - beam particles. We don't have those, yet. If they ever become part of initial
   state moduls, those module writers MUST set the status to 4, all we do here
   is respect that code (for the future).
   However, this particle won't be a parton, so we only check that for hadrons.
   In practice, this will probably require a revamp of the graph structure, both
   in the framework and in here. Then, probably also use
 GenEvent::add_beam_particle()
 - decayed particle. We don't have those either yet, but it's a feature that may
   come pretty soon.
   Here, we use that exclusively to mean decayed unstable hadrons,
   e.g. K0S -> pi pi
   In this case, the module that produced the hadron list MUST ensure to set the
   status of K0s to 2. The pions are final particles, see below.
 - Final hadrons will all be forced to status 1
 - Partons: By default, we will accept and use any code provided by the
 framework within 11<=status<=200 and use the absolute value. Parton exceptions:
            1) If the code is not HepMC-legal (|status|<11 or >200, often 0),
               we assign 12 to most and 11 to the final partons before
 hadronization (to mimic the 1, 2 scheme) 2) If there are NO hadrons in the
 event, we assume the user would like to treat the final partons (like 11 from
 above) as final particles and assign 1
 ******************************************************************************/

#include "JetScapeWriterHepMC.h"

#include <GTL/topsort.h>

#include "GTL/node.h"
#include "HardProcess.h"
#include "JetScapeLogger.h"
#include "JetScapeSignalManager.h"

using HepMC3::Units;

namespace Jetscape {

/**
 * @brief Destructor for the JetScapeWriterHepMC class.
 *
 * This destructor checks if the writer is active by calling the GetActive()
 * method. If the writer is active, it calls the Close() method to properly
 * close the writer and release any associated resources.
 */
JetScapeWriterHepMC::~JetScapeWriterHepMC() {
  if (GetActive())
    Close();
}

/**
 * @brief Writes the header information to the HepMC event file.
 *
 * This function initializes a HepMC event with the appropriate units and sets
 * various header information such as cross-section, event weight, and heavy ion
 * properties. It includes several notes regarding the units and the
 * interpretation of certain header values.
 *
 * @see
 * https://gitlab.cern.ch/hepmc/HepMC3/blob/master/include/HepMC/GenHeavyIon.h
 */
void JetScapeWriterHepMC::WriteHeaderToFile() {
  /**
   * @note Create event here - not actually writing
   * @todo TODO: GeV seems right, but I don't think we actually measure in mm
   * Should multiply all lengths by 1e-12 probably
   */
  evt = GenEvent(Units::GEV, Units::MM);

  // Expects pb, pythia delivers mb
  auto xsec = make_shared<HepMC3::GenCrossSection>();
  xsec->set_cross_section(GetHeader().GetSigmaGen() * 1e9, 0);
  xsec->set_cross_section(GetHeader().GetSigmaGen() * 1e9,
                          GetHeader().GetSigmaErr() * 1e9);
  evt.set_cross_section(xsec);
  evt.weights().push_back(GetHeader().GetEventWeight());

  auto heavyion = make_shared<HepMC3::GenHeavyIon>();
  /**
   * @note see
   * https://gitlab.cern.ch/hepmc/HepMC3/blob/master/include/HepMC/GenHeavyIon.h
   */
  if (GetHeader().GetNpart() > -1) {
    /// @note Not clear what the difference is...
    heavyion->Ncoll_hard = GetHeader().GetNcoll();
    heavyion->Ncoll = GetHeader().GetNcoll();
  }
  if (GetHeader().GetNcoll() > -1) {
    /**
     * @note Hepmc separates into target and projectile.
     * Set one? Which? Both? half to each? setting projectile for now.
     * setting both might lead to weird problems when they get added up
     */
    heavyion->Npart_proj = GetHeader().GetNpart();
  }
  if (GetHeader().GetTotalEntropy() > -1) {
    /**
     * @note nothing good in the HepMC standard. Something related to
     * mulitplicity would work
     */
  }

  if (GetHeader().GetEventPlaneAngle() > -999) {
    heavyion->event_plane_angle = GetHeader().GetEventPlaneAngle();
  }

  evt.set_heavy_ion(heavyion);

  /// @note reset per-event state
  hashadrons = false;
  hadronizationvertex = nullptr;
  pendingHadrons.clear();
  hadronsByLabel.clear();
  hadronDecayVertices.clear();
}

void JetScapeWriterHepMC::WriteEvent() {
  VERBOSE(1) << "Run JetScapeWriterHepMC: Write event # " << GetCurrentEvent();

  // create all the hepmc particles and make them available by hash for fast
  // lookup when creating the decay vertices 
  for (const auto &hadron : pendingHadrons) {

    // create the hepmc particle
    auto hepmc = castHadronToHepMC(hadron);

    // set particle hepmc status codes: 1 final state, 2 decayed, 4 special/beam
    hepmc->set_status(mapHadronStatusForHepMC(*hadron));

    // hash to hepmc hadron pointer for fast lookup
    hadronsByLabel[hadron->plabel()] = hepmc;
  }

  // create the vertices to graph mother / daughter relationships
  for (const auto &hadron : pendingHadrons) {
    auto it = hadronsByLabel.find(hadron->plabel());
    if (it == hadronsByLabel.end()) {
      continue;
    }

    // it->first is lookup key, it->second is particle pointer
    auto hepmc = it->second;

    int mother1 = hadron->mother1_label();
    int mother2 = hadron->mother2_label();

    // valid mothers are positive indexes
    if (mother1 <= 0) {
      mother1 = -1;
    }
    if (mother2 <= 0 || mother2 == mother1) {
      mother2 = -1;
    }

    // if particle has no mothers, attach to the default hadronization vertex.
    if (mother1 < 0 && mother2 < 0) {
      if (hadronizationvertex && !hepmc->production_vertex()) {
        hadronizationvertex->add_particle_out(hepmc);
      }
      continue;
    }

    int key1 = mother1;
    int key2 = mother2;

    // ensure all single mothers are in the form (key, -1)
    // and not in the form (-1, key)
    if (key1 < 0 && key2 >= 0) {
      key1 = key2;
      key2 = -1;
    }

    // ensure that for two mothers, the smaller key is always first,
    // so the same two mother indexes (regardless of order) hash to 
    // the same decay vertex
    if (key1 >= 0 && key2 >= 0 && key2 < key1) {
      std::swap(key1, key2);
    }

    // look up the decay vertex for this mother pair, if it already
    // exists, such as when adding this particle as a second daughter
    // of the same mother(s)
    auto vtxKey = std::make_pair(key1, key2);
    auto vIt = hadronDecayVertices.find(vtxKey);
    HepMC3::GenVertexPtr decayVtx = nullptr;
    if (vIt != hadronDecayVertices.end()) {
      decayVtx = vIt->second;
    }

    // lambda function to look up a mother and get its pointer
    auto findMotherParticle = [&](int label) -> HepMC3::GenParticlePtr {
      if (label < 0) {
        return nullptr;
      }
      auto mIt = hadronsByLabel.find(label);
      if (mIt == hadronsByLabel.end()) {
        return nullptr;
      }
      return mIt->second;
    };

    auto motherParticle1 = findMotherParticle(key1);
    auto motherParticle2 = findMotherParticle(key2);

    // lambda function to choose an existing decay vertex from the mothers,
    // if it exists, and warn if they have different decay vertices
    auto chooseExistingDecayVertex = [&](const HepMC3::GenParticlePtr &m,
                                         HepMC3::GenVertexPtr &selectedVtx) {

      // if mother doesn't exist or doesn't have a decay vertex, return
      if (!m || !m->end_vertex()) {
        return;
      }

      // if there's no selected decay vertex yet, use this mother's decay vertex
      if (!selectedVtx) {
        selectedVtx = m->end_vertex();
        return;
      }

      // if there's already a selected decay vertex, stick with it but warn
      // if the mothers have different daughter pairs. This might happen if
      // particle A has daughters C and D, and particle B has daughters D and E.
      // When processing daughter D, we don't know whether to attach it to 
      // the end vertex of A or B, so we just pick one (e.g. A) and warn.
      // HepMC only supports one end (decay) vertex per particle.
      if (selectedVtx != m->end_vertex()) {
        JSWARN << "Hadron " << hadron->plabel()
               << " has multiple mothers with different decay vertices.";
      }
    };

    // if no decay vertex has been created yet, create one now
    // and ad it to the vertices list. Also add it to the map.
    if (!decayVtx) {
      chooseExistingDecayVertex(motherParticle1, decayVtx);
      chooseExistingDecayVertex(motherParticle2, decayVtx);

      if (!decayVtx) {
        HepMC3::FourVector vtxPosition(hadron->x_in().x(), hadron->x_in().y(),
                                       hadron->x_in().z(), hadron->x_in().t());
        decayVtx = make_shared<GenVertex>(vtxPosition);
        vertices.push_back(decayVtx);
      }

      hadronDecayVertices[vtxKey] = decayVtx;
    }

    
    // lambda functon to attach a mother particle to the decay vertex
    // if it's not already attached.
    auto attachMotherIfCompatible = [&](const HepMC3::GenParticlePtr &mother,
                                        int label) {
      if (!mother) {
        if (label >= 0) {
          JSWARN << "Hadron " << hadron->plabel()
                 << " has a mother label " << label
                 << " but found no corresponding mother particle.";
        }
        return;
      }


      if (mother->end_vertex() && mother->end_vertex() != decayVtx) {
        JSWARN << "Mother hadron " << label
               << " has conflicting decay vertices. This hadron " << hadron->plabel()
               << " has a mother with a different decay vertex; skip attaching this mother.";
        return;
      }

      if (!mother->end_vertex()) {
        decayVtx->add_particle_in(mother);
      }
    };

    attachMotherIfCompatible(motherParticle1, key1);
    attachMotherIfCompatible(motherParticle2, key2);

    // attach this hadron as a daughter of the decay vertex
    // if it is not already attached
    if (!hepmc->production_vertex()) {
      decayVtx->add_particle_out(hepmc);
    } else if (hepmc->production_vertex() != decayVtx) {
      JSWARN << "Hadron " << hadron->plabel()
             << " already has a different production vertex; keeping existing "
                "vertex.";
    }
  }

  // Have collected all vertices now.
  // Add all vertices to the event
  for (auto v : vertices) {
    evt.add_vertex(v);
  }

  VERBOSE(1) << " found " << vertices.size() << " vertices in the list";

  // If there are no hadrons, promote final partons
  // The graph support of hepmc is a bit rudimentary.
  // easiest is to just check all childless particles
  // Note, one could just check for status==11,
  // but modules are allowed to assign that number to non-final partons
  if (!hashadrons) {
    VERBOSE(1) << " found no hadrons, promoting final partons to status 1";
    for (auto p : evt.particles()) {
      if (p->children().size() == 0) {
        if (p->status() != 11) {
          JSWARN << "Found a final parton with status!=11 : status="
                 << p->status() << ". This should not happen";
        }
        p->set_status(1);
      }
    }
  }
  evt.set_event_number(GetCurrentEvent());
  write_event(evt);
  vertices.clear();
  hadronizationvertex = 0;
  pendingHadrons.clear();
  hadronsByLabel.clear();
  hadronDecayVertices.clear();
}

/**
 * @brief This function dumps the particles in a specific parton shower to
 * the event.
 *
 * This function converts the parton shower information into a HepMC format,
 * ensuring that the topological order is respected. It uses a topological
 * sort to process the vertices and edges of the parton shower graph, creating
 * HepMC vertices and particles as needed.
 *
 * @param ps A weak pointer to the PartonShower object.
 *
 * @throws std::runtime_error If the parton shower graph is not acyclic or if
 *                            there are issues with vertex or particle creation.
 */
void JetScapeWriterHepMC::Write(weak_ptr<PartonShower> ps) {
  shared_ptr<PartonShower> pShower = ps.lock();
  if (!pShower)
    return;

  /**
   * Need topological order, see
   * https://hepmc.web.cern.ch/hepmc/differences.html
   * That means if parton p1 comes into vertex v, and p2 goes out of v,
   * then p1 has to be created (bestowed an id) before p2
   * Take inspiration from hepmc3.0.0/interfaces/pythia8/src/Pythia8ToHepMC3.cc
   * But pythia showers are different from our existing graph structure,
   * So instead try to modify the first attempt to respect top. order
   * and don't create vertices and particles more than once
   *
   * Using GTL's topsort
   * 1. Check that our graph is sane
   */
  if (!pShower->is_acyclic())
    throw std::runtime_error(
        "PROBLEM in JetScapeWriterHepMC: Graph is not acyclic.");

  /// @note 2.
  topsort topsortsearch;
  topsortsearch.scan_whole_graph(true);
  topsortsearch.start_node(); /** @note defaults to first node */
  topsortsearch.run(*pShower);
  /// @note this is a topsort::topsort_iterator
  auto nEnd = topsortsearch.top_order_end();

  /// @note Need to keep track of already created ones
  map<int, GenParticlePtr> CreatedPartons;

  /// @note probably unnecessary, used for consistency checks
  map<int, bool> vused;
  bool foundRoot = false;
  for (auto nIt = topsortsearch.top_order_begin(); nIt != nEnd; ++nIt) {
    // cout << *nIt << "  " << nIt->indeg() << "  " << nIt->outdeg() << endl;

    /// @note Should be the only time we see this node.
    if (vused.find(nIt->id()) != vused.end()) {
      throw std::runtime_error(
          "PROBLEM in JetScapeWriterHepMC: Reusing a vertex.");
    }
    vused[nIt->id()] = true;

    /**
     * @note 0. No incoming edges?
     * ---------------------
     * This is typically a shower initiator.
     * HepMC needs an incoming and outgoing particle for every vertex.
     * That could be a place to attach partons or ions.
     * We could also do both, but as of now, JETSCAPE actually attaches a dummy
     * vertex to the start of the initiator, that can safely go away.
     * Previously, we attached a dummy or clone of the outgoing one here.
     * Instead, we can just skip the vertex. Its outgoing edges will be picked
     * up as incomers in a later vertex. Note that the [0]=>[1] connection in
     * JETSCAPE already uses a dummy node[0], and [1] is at time t=0; removing
     * that seems correct.
     */
    if (nIt->indeg() == 0)
      continue;

    /**
     * @note 1. Create a new vertex.
     * --------------------------------------------
     */
    auto v = castVtxToHepMC(pShower->GetVertex(*nIt));

    /**
     * @note 2. Incoming edges
     * --------------------------------------------
     * In the current framework, it should only be one.
     * So we will catch anything more but provide a mechanism that should
     * work anyway.
     */
    if (nIt->indeg() > 1) {
      JSWARN << "Found more than one mother parton! Should only happen if we "
                "added medium particles. "
             << "The code should work, but proceed with caution";
    }

    auto inIt = nIt->in_edges_begin();
    auto inEnd = nIt->in_edges_end();
    for (/* nop */; inIt != inEnd; ++inIt) {
      auto phepin = CreatedPartons.find(inIt->id());
      if (phepin != CreatedPartons.end()) {
        /// @note We should already have one!
        v->add_particle_in(phepin->second);
      } else {
        // This indicates we skipped an earlier vertex without incomers.
        // JSWARN << "Incoming particle out of nowhere. This could maybe happen
        // "
        //           "if we pick up medium particles "
        //        << " but is probably a topsort problem. Try using the code "
        //           "after this throw() but be very careful.";
        // throw std::runtime_error("PROBLEM in JetScapeWriterHepMC: Incoming "
        //                          "particle out of nowhere.");

        auto in = pShower->GetParton(*inIt);
        auto hepin = castPartonToHepMC(in);
        auto status = std::abs(hepin->status());
        if (status < 11 || status > 200) {
          /// @note incoming edge can't be final
          status = 12;
        }
        hepin->set_status(status);
        CreatedPartons[inIt->id()] = hepin;
        v->add_particle_in(hepin);

        if (nIt->outdeg() == 0) {
          /**
           * @note However, motherless AND childless particles do exist
           * I.e., a shower initiator that never actually showers
           * For this, we need an out going clone, much like 3) below
           */
          auto hepout = castPartonToHepMC(in);
          /**
           * @note Since the status information is preserved in the incomer
           * we'll force 11 Note: if we later see in WriteEvent() that there
           * are no hadrons, this will be overwritten to 1
           */
          hepout->set_status(11);
          /**
           * @note Note that we do not register this particle. Since it's
           * pointing nowhere, it can never be reused.
           */
          v->add_particle_out(hepout);
        }
      }
    }

    /**
     * @note 3. Outgoing edges?
     * --------------------------------------------
     * 3.1: No. Need to create one.
     * We'll use this opportunity to copy the incomer but give it a final code
     */
    if (nIt->outdeg() == 0) {
      if (nIt->indeg() != 1) {
        /**
         * @note This won't work with multiple incomers
         * (but that's pretty unphysical)
         */
        throw std::runtime_error(
            "PROBLEM in JetScapeWriterHepMC: Need exactly "
            "one parent to clone final state partons.");
      }
      auto in = pShower->GetParton(*(nIt->in_edges_begin()));
      auto hepout = castPartonToHepMC(in);
      /**
       * @note an outgoing edge without terminator is "final"
       * Since the status information is preserved in the incomer, we'll force
       * 11 Note: if we later see in WriteEvent() that there are no hadrons,
       * this will be overwritten to 1
       */
      hepout->set_status(11);
      v->add_particle_out(hepout);
      /**
       * @note Note that we do not register this particle. Since it's pointing
       * nowhere, it can never be reused.
       */
    }

    /// @note 3.2: Otherwise use and register the outgoing edge
    if (nIt->outdeg() > 0) {
      auto outIt = nIt->out_edges_begin();
      auto outEnd = nIt->out_edges_end();
      for (/* nop */; outIt != outEnd; ++outIt) {
        if (CreatedPartons.find(outIt->id()) != CreatedPartons.end()) {
          throw std::runtime_error(
              "PROBLEM in JetScapeWriterHepMC: Trying to "
              "recreate a preexisting GenParticle.");
        }
        auto out = pShower->GetParton(*outIt);
        auto hepout = castPartonToHepMC(out);
        if (!hepout->status()) {
          /// @note incoming and outgoing -> status 12
          hepout->set_status(12);
        }

        CreatedPartons[outIt->id()] = hepout;
        v->add_particle_out(hepout);
      }
    }

    vertices.push_back(v);
  }
}

/**
 * @brief Writes a hadron to the HepMC3 format.
 *
 * This function takes a weak pointer to a Hadron object, locks it to obtain a
 * shared pointer, and writes it to the HepMC3 format. If the hadronization
 * vertex does not exist, it creates one with a dummy position and a dummy
 * mother particle. All hadrons are attached to this hadronization vertex. If
 * the hadron's status is not specified, it is set to 1 by default.
 *
 * @param h A weak pointer to a Hadron object.
 */
void JetScapeWriterHepMC::Write(weak_ptr<Hadron> h) {
  auto hadron = h.lock();
  if (!hadron)
    return;

  /**
   * @note No clear source for most hadrons
   * Also, a graph with e.g. recombination hadrons would have loops,
   * (though the direction should still make it acyclic?)
   * Not sure how this is supposed to be done in HepMC3
   * Our solution: Attach all hadrons to one dedicated hadronization vertex.
   * Future option: Have separate shower and bulk vertices?
   */

  /// @note Create if it doesn't exist yet
  if (!hadronizationvertex) {
    /// @note dummy position, set it to a late time...
    HepMC3::FourVector vtxPosition(0, 0, 0, 100);
    hadronizationvertex = make_shared<GenVertex>(vtxPosition);

    /**
     * @note dummy mother -- could also maybe use the first/hardest shower
     * initiator
     */
    HepMC3::FourVector pmom(0, 0, 0, 0);
    make_shared<GenParticle>(pmom, 0, 0);
    hadronizationvertex->add_particle_in(make_shared<GenParticle>(pmom, 0, 0));

    vertices.push_back(hadronizationvertex);
    hashadrons = true;
  }

  pendingHadrons.push_back(hadron);
}

/**
 * @brief Initializes the JetScape HepMC Writer.
 *
 * This function checks if the writer is active. If it is active, it logs
 * an informational message indicating that the JetScape HepMC Writer has
 * been initialized along with the output file name.
 */
void JetScapeWriterHepMC::Init() {
  if (GetActive()) {
    JSINFO << "JetScape HepMC Writer initialized with output file = "
           << GetOutputFileName();
  }
}

/**
 * @brief Executes the JetScapeWriterHepMC process.
 */
void JetScapeWriterHepMC::Exec() {
  /**
   * @note This function currently does not perform any operations.
   */
}
}  // end namespace Jetscape
