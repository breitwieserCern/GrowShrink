// -----------------------------------------------------------------------------
//
// Copyright (C) The BioDynaMo Project.
// All Rights Reserved.
//
// Licensed under the Apache License, Version 2.0 (the "License");
// you may not use this file except in compliance with the License.
//
// See the LICENSE file distributed with this work for details.
// See the NOTICE file distributed with this work for additional information
// regarding copyright ownership.
//
// -----------------------------------------------------------------------------
#ifndef SHRINKGROW2_H_
#define SHRINKGROW2_H_

#include "biodynamo.h"
//#include "diffusion_biology_modules.h"
//#include "hello2.cc"


namespace bdm {

// 0. extending cell behaviour "class"
BDM_SIM_OBJECT(MyCell, Cell){

  BDM_SIM_OBJECT_HEADER(MyCell,Cell,1,cell_VC_,cell_IR_,cell_CR_,cell_T_,cell_L_,cell_IL_);
   public:
    MyCellExt() {}
    explicit MyCellExt(const std::array<double, 3>& position) : Base(position) {}
// Setter for dV/dt
    void SetGrowthRate( double L,double A,double R,double T,double O) {
      cell_VC_[kIdx] = L*A*R*T*O;
    }
    double GetGrowthRate() const { return cell_VC_[kIdx]; }
// Setter for Initial Radius of cell
    void SetIR( double IR ) {
      cell_IR_[kIdx] = IR;
    }
    double GetIR() const { return cell_IR_[kIdx]; }
// Setter for Cooling Rate
    void SetCR( double CR ) {
      cell_CR_[kIdx] = CR;
    }
    double GetCR() const { return cell_CR_[kIdx]; }
// Temperature Setter
   void SetT( double T ) {
      cell_T_[kIdx] = T;
    }
    double GetT() const { return cell_T_[kIdx]; }
// Hydraulic permeability setter for T = const , this is also constant.
   void SetL(double LR, double AE ,double R ,double TI, double RT){
      cell_L_[kIdx] = LR*exp((-AE/R)*((1/TI)-(1/RT)));

  }
      double GetL() const { return cell_T_[kIdx]; }
//Initial Lp setter
    void SetIL( double IL ) {
      cell_IL_[kIdx] = IL;
    }
    double GetIL() const { return cell_IL_[kIdx]; }


   private:
    vec<int> cell_VC_;
    vec<int> cell_IR_;
    vec<int> cell_CR_;
    vec<int> cell_T_;
    vec<int> cell_L_;
    vec<int> cell_IL_;
};

// 1. Growth behaviour
struct GrowthModule : public BaseBiologyModule {
  GrowthModule() : BaseBiologyModule(gAllEventIds) {}
  GrowthModule(double threshold, double growth_rate,
             std::initializer_list<EventId> event_list)
      : BaseBiologyModule(event_list),
        threshold_(threshold),
        growth_rate_(growth_rate) {}

  /// Default event constructor
  template <typename TEvent, typename TBm>
  GrowthModule(const TEvent& event, TBm* other, uint64_t new_oid = 0) {
    threshold_ = other->threshold_;
    growth_rate_ = other->growth_rate_;
  }


  template <typename T, typename TSimulation = Simulation<>>
  void Run(T* cell) {

    if (cell->GetGrowthRate() > 0) {
      cell->ChangeVolume(cell->GetGrowthRate()*1); // Decreasing volume

    } else if (cell->GetGrowthRate() < 0) {
      double IR = cell->GetIR();
      double MT = 263;
        if (cell->GetDiameter() > 0.7*IR){
          cell->ChangeVolume(cell->GetGrowthRate()*1); //Increasing volume
          if (cell -> GetT() > MT){
             double CT = cell -> GetT();
             cell->SetT(CT - cell->GetCR());
             double IL = cell->GetIL();
             cell->SetL(IL,1,1,CT,274);
          } else {
            cell->SetT(cell->GetT());
          }
        } else {
          cell->ChangeVolume(1); // Cell has reached Osmal equilibrium
        }
    } else {
      cell->ChangeVolume(1); // Osmoality is balanced therefoe volume is constant.
    }
  }

 private :
  BDM_CLASS_DEF_NV(GrowthModule, 1);
  double threshold_ = 0;
  double growth_rate_  = 1; // NEED TO CORRECT;
};


// 2 . Define compile time parameter
BDM_CTPARAM() {
  BDM_CTPARAM_HEADER();
  using SimObjectTypes = CTList<MyCell>;  // use MyCell object

  // Override default BiologyModules for Cell
  BDM_CTPARAM_FOR(bdm, MyCell) { using BiologyModules = CTList<GrowthModule>; };
};

inline int Simulate(int argc, const char** argv) {
  auto set_param = [](auto* param) {
    param->bound_space_ = true;
    param->min_bound_ = 0;
    param->max_bound_ = 1000;  // cube of 100*100*100
    param->run_mechanical_interactions_ = true;
  };

  Simulation<> simulation(argc, argv, set_param);
  auto* rm = simulation.GetResourceManager();
  auto* param = simulation.GetParam();
  auto* random = simulation.GetRandom();

  size_t nb_of_cells = 1000;  // number of cells in the simulation
  double x_coord, y_coord, z_coord;

  // create a structure to contain cells
  auto* cells = rm->template Get<MyCell>();
  // allocate the correct number of cell in our cells structure before
  // cell creation
  cells->reserve(nb_of_cells);

  for (size_t i = 0; i < nb_of_cells; ++i) {
    // our modelling will be a cell cube of 100*100*100
    // random double between 0 and 100
    x_coord = random->Uniform(param->min_bound_, param->max_bound_);
    y_coord = random->Uniform(param->min_bound_, param->max_bound_);
    z_coord = random->Uniform(param->min_bound_, 500);

    // creating the cell at position x, y, z
    MyCell cell({x_coord, y_coord, z_coord});
    // set cell parameters
    cell.SetDiameter(40);
    cell.SetIR(cell->GetDiameter());
    cell.SetCR(1);
    cell.SetT(274);
    cell.SetIL(20);
    cell.SetL(20,1,1,274,274);
    cell.SetGrowthRate(cell->GetL(),1,1,cell->GetT() - cell->GetCR(),-1);
    cell.AddBiologyModule(GrowthModule());
    cells->push_back(cell);  // put the created cell in our cells structure
  }

    for (size_t i = 0; i < nb_of_cells; ++i) {
    // our modelling will be a cell cube of 100*100*100
    // random double between 0 and 100
    x_coord = random->Uniform(param->min_bound_, param->max_bound_);
    y_coord = random->Uniform(param->min_bound_, param->max_bound_);
    z_coord = random->Uniform(500, param->max_bound_);

    // creating the cell at position x, y, z
    MyCell cell({x_coord, y_coord, z_coord});
    cell.SetDiameter(30);
    cell.SetIR(cell->GetDiameter());
    cell.SetCR(1);
    cell.SetT(274);
    cell.SetIL(30);
    cell.SetL(30,1,1,274,274);
    cell.SetGrowthRate(cell->GetL(),1,1,cell->GetT() - cell->GetCR(),-1);
    cell.AddBiologyModule(GrowthModule());
    cells->push_back(cell);  // put the created cell in our cells structure
    }
      for (size_t i = 0; i < nb_of_cells; ++i) {
    // our modelling will be a cell cube of 100*100*100
    // random double between 0 and 100
    x_coord = random->Uniform(param->min_bound_, param->max_bound_);
    y_coord = random->Uniform(param->min_bound_, param->max_bound_);
    z_coord = random->Uniform(250, 750);

    // creating the cell at position x, y, z
    MyCell cell({x_coord, y_coord, z_coord});
    cell.SetDiameter(35);
    cell.SetAdherence(30);
    cell.SetIR(cell->GetDiameter());
    cell.SetCR(1);
    cell.SetT(274);
    cell.SetIL(25);
    cell.SetL(25,1,1,274,274);
    cell.SetGrowthRate(cell->GetL(),1,1,cell->GetT() - cell->GetCR(),-1);
    cell.AddBiologyModule(GrowthModule());
    cells->push_back(cell);  // put the created cell in our cells structure
    }

    cells->Commit();  // commit cells

  // Run simulation for one timestep
  simulation.GetScheduler()->Simulate(50);

  std::cout << "Simulation completed successfully! FTH" << std::endl;
  return 0;
}

}  // namespace bdm

#endif  // SHRINKGROW2_H_
