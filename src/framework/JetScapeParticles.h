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

/**
 * @brief Base particle class in JETSCAPE, derived privately from FastJet
 * PseudoJet.
 *
 * Provides `JetScapeParticleBase` and derived classes (`Parton`, `Photon`,
 * `Hadron`).
 *
 * @class Jetscape::JetScapeParticleBase
 *
 * A `JetScapeParticleBase` derives privately from FastJet `PseudoJet` and
 * adds additional information:
 * - PID (from PDG) and rest mass (should eventually be coupled, keeping only
 * PID internally)
 * - A location (creation point) 4-vector
 * - A label and a status code
 * - Additional information (some of which should be moved to derived classes or
 * `UserInfo`)
 *
 * ### Design choice
 * Protected inheritance is used due to mismatches between conventions:
 * - Most of the theory community expects the 0-component to be time/energy.
 * - FastJet (and packages like ROOT) instead use the 4th component for
 * time/energy.
 *
 * Private inheritance allows inheriting safe methods (via `using` in C++11),
 * while protecting users from unsafe component access.
 * This is necessary to prevent unwanted use of constructors and getters that
 * assume FastJet’s indexing conventions.
 * If these were removed or adapted, we could use public inheritance and obtain
 * a proper "is-a" relationship.
 *
 * ### Usage
 * The base class can be used directly, but it is recommended to use derived
 * classes such as:
 * - `Parton`
 * - `Hadron`
 * - `Photon`
 * - (future) `Lepton`
 *
 * @warning
 * `PseudoJet` has no concept of rest mass. Its mass-related functions assume:
 * \f[
 * M^2 = E^2 - \vec{p}^2
 * \f]
 * For off-shell partons, the correct interpretation is:
 * \f[
 * M^2 = E^2 - \vec{p}^2 = M_0^2 + Q^2 = M_0^2 + t
 * \f]
 * This creates the risk that functions like `mass()` or `reset_PtYPhiM(...)`
 * will give misleading results. To protect against this, "mass"-related
 * functions in `PseudoJet` are intentionally not made available.
 *
 * ### Future considerations
 * - Making a Pythia8 installation mandatory: with guaranteed Pythia, rest mass
 *   lookup could use PDG data automatically.
 * - If ROOT were mandatory, `TLorentzVector` could replace the custom
 * `FourVector`.
 * - If HepMC were mandatory, `HepMC::FourVector` could replace the custom
 * `FourVector`.
 */
#ifndef JETSCAPEPARTICLES_H
#define JETSCAPEPARTICLES_H

// clang-format tool is disabled for these includes because clang-format
// will alphabetize them, but their included order matters
// clang-format off
#include <stdio.h>
#include <math.h>
#include "GTL/graph.h"
#include "JetScapeConstants.h"
#include "FourVector.h"
#include "fjcore.hh"
#include "JetScapeLogger.h"
#include "PartonShower.h"

#include "Pythia8/Pythia.h"

#include <vector>
#include <memory>
#include <iostream>
#include <sstream>
#include <iomanip>
// clang-format on

using std::ostream;
using std::weak_ptr;

namespace Jetscape {

class PartonShower;

/**
 * @class JetScapeParticleBase
 * @brief Base particle class in JETSCAPE, deriving from FastJet's PseudoJet.
 *
 * Adds physics-specific information to `fjcore::PseudoJet`, including:
 * - PDG ID and rest mass
 * - Status and label
 * - Production vertex (`x_in_`)
 * - An optional jet velocity vector (`jet_v_`)
 * - Control information for responsible modules (e.g. energy loss)
 *
 * Private inheritance is used to restrict unsafe direct access to
 * FastJet components (notably time vs. energy convention mismatches).
 *
 * @warning Mass-related functions from `PseudoJet` are deliberately
 * hidden to avoid misuse in off-shell contexts.
 */
class JetScapeParticleBase : protected fjcore::PseudoJet {
  friend class fjcore::PseudoJet;
  // unsafe
  // using fjcore::PseudoJet::PseudoJet;
  // using fjcore::PseudoJet::operator() (int i) const ;
  // inline double operator [] (int i) const { return (*this)(i); }; // this too

 public:
  // Disallow reset, it assumes different logic
  // using fjcore::PseudoJet::reset(double px, double py, double pz, double E);
  // using fjcore::PseudoJet::reset(const PseudoJet & psjet) {
  //   inline void reset_momentum(double px, double py, double pz, double E);
  //   inline void reset_momentum(const PseudoJet & pj);

  inline void reset_momentum(const double px, const double py, const double pz,
                             const double e) {
    fjcore::PseudoJet::reset_momentum(px, py, pz, e);
  }

  inline void reset_momentum(const FourVector &p) {
    fjcore::PseudoJet::reset_momentum(p.x(), p.y(), p.z(), p.t());
  }

  // Disallow the valarray return.
  // Can replace/and or provide a double[4] version with the right assumptions
  // std::valarray<double> four_mom() const;
  // enum { X=0, Y=1, Z=2, T=3, NUM_COORDINATES=4, SIZE=NUM_COORDINATES };

  // This _should_ work, but not taking chances for now
  // PseudoJet & boost(const PseudoJet & prest);
  // PseudoJet & unboost(const PseudoJet & prest);

  // Replace with appropriate functions. m is a tricky one.
  // inline double  m2() const {return (_E+_pz)*(_E-_pz)-_kt2;}
  // inline double  m() const;
  // inline double mperp2() const {return (_E+_pz)*(_E-_pz);}
  // inline double mperp() const {return sqrt(std::abs(mperp2()));}
  // inline double mt2() const {return (_E+_pz)*(_E-_pz);}
  // inline double mt() const {return sqrt(std::abs(mperp2()));}

  // Disallow functions containing M
  // inline void reset_PtYPhiM(...);
  // void reset_momentum_PtYPhiM(double pt, double y, double phi, double m=0.0);

  // void set_cached_rap_phi(double rap, double phi);

  /// No implicit cast to PseudoJet is allowed, provide a conversion
  fjcore::PseudoJet GetPseudoJet() const { return PseudoJet(*this); }

  // import safe functions
  using fjcore::PseudoJet::e;
  using fjcore::PseudoJet::E;
  using fjcore::PseudoJet::px;
  using fjcore::PseudoJet::py;
  using fjcore::PseudoJet::pz;

  using fjcore::PseudoJet::eta;
  using fjcore::PseudoJet::kt2;
  using fjcore::PseudoJet::perp;
  using fjcore::PseudoJet::perp2;
  using fjcore::PseudoJet::phi;
  using fjcore::PseudoJet::phi_02pi;
  using fjcore::PseudoJet::phi_std;
  using fjcore::PseudoJet::pseudorapidity;
  using fjcore::PseudoJet::pt;
  using fjcore::PseudoJet::pt2;
  using fjcore::PseudoJet::rap;
  using fjcore::PseudoJet::rapidity;

  using fjcore::PseudoJet::Et;
  using fjcore::PseudoJet::Et2;
  using fjcore::PseudoJet::modp;
  using fjcore::PseudoJet::modp2;

  using fjcore::PseudoJet::beam_distance;
  using fjcore::PseudoJet::delta_phi_to;
  using fjcore::PseudoJet::delta_R;
  using fjcore::PseudoJet::kt_distance;
  using fjcore::PseudoJet::plain_distance;
  using fjcore::PseudoJet::squared_distance;

  using fjcore::PseudoJet::operator*=;
  using fjcore::PseudoJet::operator/=;
  using fjcore::PseudoJet::operator+=;
  using fjcore::PseudoJet::operator-=;

  using fjcore::PseudoJet::InexistentUserInfo;
  using fjcore::PseudoJet::set_user_index;
  using fjcore::PseudoJet::user_index;
  using fjcore::PseudoJet::UserInfoBase;

  using fjcore::PseudoJet::has_user_info;
  using fjcore::PseudoJet::set_user_info;
  using fjcore::PseudoJet::user_info;
  using fjcore::PseudoJet::user_info_ptr;
  using fjcore::PseudoJet::user_info_shared_ptr;

  using fjcore::PseudoJet::description;
  // In principle, these might be okay, but ClusterSequences should
  // be made after explicitly transforming to a proper PseudoJet
  // using fjcore::PseudoJet::has_associated_cluster_sequence;
  // using fjcore::PseudoJet::has_associated_cs;
  // using fjcore::PseudoJet::has_valid_cluster_sequence;
  // using fjcore::PseudoJet::has_valid_cs;
  // using fjcore::PseudoJet::associated_cluster_sequence;
  // using fjcore::PseudoJet::associated_cs;
  // using fjcore::PseudoJet::validated_cluster_sequence;
  // using fjcore::PseudoJet::validated_cs;
  // using fjcore::PseudoJet::set_structure_shared_ptr;
  // using fjcore::PseudoJet::has_structure;
  // using fjcore::PseudoJet::structure_ptr;
  // using fjcore::PseudoJet::structure_non_const_ptr;
  // using fjcore::PseudoJet::validated_structure_ptr;
  // using fjcore::PseudoJet::structure_shared_ptr;
  // ... more

 public:
  /// Default constructor
  JetScapeParticleBase() : PseudoJet(){};

  /// Construct particle from label, ID, status, momentum and position
  JetScapeParticleBase(int label, int id, int stat, const FourVector &p,
                       const FourVector &x);
  /// Construct from label, ID, status, kinematics and optional position
  JetScapeParticleBase(int label, int id, int stat, double pt, double eta,
                       double phi, double e, double *x = 0);
  /// Construct with explicit rest mass
  JetScapeParticleBase(int label, int id, int stat, const FourVector &p,
                       const FourVector &x, double mass);
  /// Copy constructor
  JetScapeParticleBase(const JetScapeParticleBase &srp);

  /// Destructor
  virtual ~JetScapeParticleBase();

  /// Reset internal state
  void clear();

  // --- Setters ---
  void set_label(int label);     ///< Set label (event record line)
  void set_id(int id);           ///< Set PDG ID
  void set_stat(int stat);       ///< Set status code
  void set_x(double x[4]);       ///< Set production position
  void init_jet_v();             ///< Initialize jet_v
  void set_jet_v(double v[4]);   ///< Set jet velocity vector
  void set_jet_v(FourVector j);  ///< Set jet velocity vector

  /** Set responsible module (e.g. energy loss).
   * @return false if already controlled. */
  bool SetController(string controller = "") {
    bool wascontrolled = controlled_;
    controlled_ = true;
    controller_ = controller;
    return wascontrolled;
  };
  /// Relinquish responsibility
  void UnsetController() {
    controller_ = "";
    controlled_ = false;
  };

  // --- Getters ---
  const int pid() const;      ///< Get PDG ID
  const int pstat() const;    ///< Get status code
  const int plabel() const;   ///< Get label
  const double time() const;  ///< Get time component of position
  std::vector<JetScapeParticleBase> parents();  ///< Retrieve parents
  const FourVector p_in() const;                ///< Get incoming momentum
  const FourVector &x_in() const;               ///< Get production position
  const FourVector &jet_v() const;              ///< Get jet velocity vector

  const double restmass();  ///< Get rest mass
  const double p(int i);    ///< Get momentum component
  double pl();
  const double nu();
  const double t_max();

  /// Assignment operators
  virtual JetScapeParticleBase &operator=(JetScapeParticleBase &c);
  virtual JetScapeParticleBase &operator=(const JetScapeParticleBase &c);

  string GetController() const {
    return controller_;
  };  ///< Name of controlling module
  bool GetControlled() const {
    return controlled_;
  };  ///< Whether controlled by module

  /// Static Pythia instance for PDG lookup
  static Pythia8::Pythia InternalHelperPythia;

 protected:
  void set_restmass(double mass_input);  ///< Set rest mass internally

  int pid_;           ///< PDG ID
  int pstat_;         ///< Status code
  int plabel_;        ///< Label in event record
  double mass_;       ///< Rest mass of particle
  FourVector x_in_;   ///< Production vertex
  FourVector jet_v_;  ///< Jet four vector, without gamma factor (so not really
                      ///< a four vector)
  bool controlled_ = false;  ///< Whether controlled by module
  string controller_ = "";   ///< Name of controlling module
};
// END BASE CLASS

/// Stream output operator
ostream &operator<<(ostream &output, JetScapeParticleBase &p);

/**
 * @class Parton
 * @brief Derived class for partons.
 *
 * Adds parton-specific properties such as:
 * - Formation times
 * - Color/anticolor
 * - Shower membership
 */
class Parton : public JetScapeParticleBase {
 public:
  Parton(int label, int id, int stat, const FourVector &p, const FourVector &x);
  Parton(int label, int id, int stat, double pt, double eta, double phi,
         double e, double *x = 0);
  Parton(const Parton &srp);

  Parton &operator=(Parton &c);
  Parton &operator=(const Parton &c);

  // --- Physics-specific ---
  virtual void set_mean_form_time();             ///< Set mean formation time
  virtual void set_form_time(double form_time);  ///< Set formation time
  virtual double form_time();                    ///< Get formation time
  virtual const double mean_form_time();         ///< Get mean formation time
  virtual void reset_p(double px, double py,
                       double pz);  ///< Reset spatial momentum

  // --- Color ---
  virtual void set_color(unsigned int col);        ///< Set color
  virtual void set_anti_color(unsigned int acol);  ///< Set anti-color
  virtual void set_max_color(unsigned int col);    ///< Set maximum color
  virtual void set_min_color(unsigned int col);    ///< Set minimum color
  virtual void set_min_anti_color(
      unsigned int acol);         ///< Set minimum anti-color
  unsigned int color();           ///< Get color
  unsigned int anti_color();      ///< Get anticolor
  unsigned int max_color();       ///< Get max color
  unsigned int min_color();       ///< Get min color
  unsigned int min_anti_color();  ///< Get min anticolor

  bool isPhoton(int pid);  ///< Check if particle is photon

  const double t();      ///< Get virtuality
  void set_t(double t);  ///< Set virtuality of particle (rescales the spatial
                         ///< component)
  const int edgeid() const;       ///< Position in shower graph
  void set_edgeid(const int id);  ///< Set the edge id

  void set_shower(const shared_ptr<PartonShower> pShower);  ///< Set shower
  void set_shower(const weak_ptr<PartonShower> pShower);  ///< Set shower (weak)
  const weak_ptr<PartonShower> shower() const;            ///< Get shower
  std::vector<Parton> parents();                          ///< Get parents

 protected:
  double mean_form_time_;      ///< Mean formation time
  double form_time_;           ///< event by event formation time
  unsigned int Color_;         ///< Large Nc color of parton
  unsigned int antiColor_;     ///< Large Nc anti-color of parton
  unsigned int MaxColor_;      ///< the running maximum color
  unsigned int MinColor_;      ///< color of the parent
  unsigned int MinAntiColor_;  ///< anti-color of the parent

  weak_ptr<PartonShower> pShower_;  ///< shower that this parton belongs to
  int edgeid_;                      ///< Position in the shower graph

  void initialize_form_time();      ///< Initialize formation time
  void CheckAcceptability(int id);  ///< Check if parton ID is acceptable
};

/**
 * @class Hadron
 * @brief Derived class for hadrons.
 *
 * Adds hadron-specific attributes such as decay width.
 */
class Hadron : public JetScapeParticleBase {
 public:
  Hadron(int label, int id, int stat, const FourVector &p, const FourVector &x);
  Hadron(int label, int id, int stat, double pt, double eta, double phi,
         double e, double *x = 0);
  Hadron(int label, int id, int stat, const FourVector &p, const FourVector &x,
         double mass);
  Hadron(const Hadron &srh);

  Hadron &operator=(Hadron &c);
  Hadron &operator=(const Hadron &c);

  void set_mother_labels(int mother1, int mother2) {
    mother1_label_ = mother1;
    mother2_label_ = mother2;
  }
  void set_daughter_labels(int daughter1, int daughter2) {
    daughter1_label_ = daughter1;
    daughter2_label_ = daughter2;
  }
  
  // Getters for mother and daughter labels
  int mother1_label() const { return mother1_label_; }
  int mother2_label() const { return mother2_label_; }
  int daughter1_label() const { return daughter1_label_; }
  int daughter2_label() const { return daughter2_label_; }

  void set_decay_width(double width) { width_ = width; }  ///< Set decay width
  double decay_width() { return (width_); }               ///< Get decay width

  /// Hadron may be used to handle electrons, gammas, ... as well
  /// In addition, not all generated ids may be in the database
  /// Currently, we add these manually. Could also reject outright.
  bool CheckOrForceHadron(const int id, const double mass = 0);

  /// Check if hadron has no spatial position (all components zero)
  bool has_no_position();

 protected:
  double width_;  ///< Decay width
  int mother1_label_ = -1; ///< Label of first mother particle
  int mother2_label_ = -1; ///< Label of second mother particle
  int daughter1_label_ = -1; ///< Label of first daughter particle
  int daughter2_label_ = -1; ///< Label of second daughter particle
};

/**
 * @class Photon
 * @brief Derived class for photons.
 *
 * Specialization of Parton class.
 */
class Photon : public Parton {
 public:
  Photon(int label, int id, int stat, const FourVector &p, const FourVector &x);
  Photon(int label, int id, int stat, double pt, double eta, double phi,
         double e, double *x = 0);
  Photon(const Photon &srh);

  Photon &operator=(Photon &ph);
  Photon &operator=(const Photon &ph);
};

/**
 * @class Qvector
 * @brief Flow Q-vector container.
 *
 * Stores multi-dimensional histograms of hadron distributions in
 * transverse momentum and rapidity, with configurable harmonic order.
 */
class Qvector {
 public:
  Qvector(double pt_min, double pt_max, int npt, double y_min, double y_max,
          int ny, int norder, int pid, int rapidity_type);

  void fill(double pt_in, double y_in, int col_in, double val);  ///< Fill bin
  void fill_particle(const shared_ptr<Hadron> &hadron);  ///< Fill from hadron

  int get_pdgcode() const { return pid_; }    ///< Get PDG code
  int get_npt() const { return npt_; }        ///< Get number of pT bins
  int get_ny() const { return ny_; }          ///< Get number of rapidity bins
  int get_norder() const { return norder_; }  ///< Get harmonic order
  int get_ncols() const { return ncols_; }    ///< Get number of columns
  int get_total_num() { return total_num_; }  ///< Get total entries
  double get_value(int i, int j, int k) const {
    return hist_[i][j][k];
  }  ///< Access bin content

  double get_dpt() const { return dpt_; }             ///< Bin width in pT
  double get_dy() const { return dy_; }               ///< Bin width in rapidity
  double get_pt(int idx) const;                       ///< Center of pT bin
  double get_y(int idx) const;                        ///< Center of y bin
  void set_header(std::string a);                     ///< Set output header
  std::string get_header() const { return header_; }  ///< Get output header

 private:
  double pt_min_, pt_max_, y_min_, y_max_;
  int npt_, ny_, ncols_, norder_;
  int pid_, rapidity_type_, total_num_;
  double dpt_, dy_;
  std::vector<std::vector<std::vector<double>>> hist_;
  std::string header_;
  std::vector<double> gridpT_;
  std::vector<double> gridy_;
};

};  // namespace Jetscape

#endif  // JETSCAPEPARTICLES_H
