//
//  AO_sim.cpp
//  K_Wave_C
//
//  Created by jacob on 12/3/12.
//  Copyright (c) 2012 BMPI. All rights reserved.
//

#include <exception>

#include "AO_sim.h"





AO_Sim::AO_Sim()
{
    // Responsible for running the monte-carlo simulation.
    da_boost = new MC_Boost();
    
    // Responsible for running the k-Wave simulation.
    KSpaceSolver = new TKSpaceFirstOrder3DSolver();
    
    /// Set default values for private members.
    m_Laser_injection_coords.x = -1;
    m_Laser_injection_coords.y = -1;
    m_Laser_injection_coords.z = -1;
    
}



AO_Sim::~AO_Sim()
{
    if (da_boost)
        delete da_boost;
    
    if (KSpaceSolver)
        delete KSpaceSolver;
    
}



/// Run the actual acousto-optic simulation.
void
AO_Sim::Run_acousto_optics_sim(TParameters * Parameters)
{
    
    
    /// Load data from disk
    cout << "Loading k-Wave data ........";
    try {
        KSpaceSolver->LoadInputData();
    }catch (ios::failure e) {
        cout << "Failed!\nK-Wave panic: Data loading was not successful!\n%s\n"
             << e.what();
        cerr << "K-Wave panic: Data loading was not successful! \n%s\n"
             << e.what();
        exit(EXIT_FAILURE);
    }
    cout << ".... done\n";
    
    //fprintf(stdout,"Elapsed time:          %8.2fs\n",KSpaceSolver.GetDataLoadTime());
    
    /// start computation of k-Wave simulation.
    cout << ".......... k-Wave Computation ...........\n";
    KSpaceSolver->PreCompute();
    
    
//    ActPercent = 0;
    KSpaceSolver->FromMain_PrintOutputHeader();
//
    KSpaceSolver->StartIterationTime();
    for (KSpaceSolver->SetTimeIndex(0); KSpaceSolver->GetTimeIndex() < Parameters->Get_Nt(); KSpaceSolver->IncrementTimeIndex()){
        //    for (t_index = 0; t_index < Parameters->Get_Nt(); t_index ++){
        //
        //
        KSpaceSolver->FromMain_Compute_uxyz();
//        //        Compute_uxyz();
//        //
//        //        // add in the velocity u source term
        KSpaceSolver->FromMain_Add_u_source();
//        //        Add_u_source();
//        //
//        //
//        //        // add in the transducer source term (t = t1) to ux
//        //        if (Parameters->Get_transducer_source_flag() > t_index)
        if (Parameters->Get_transducer_source_flag() > KSpaceSolver->GetTimeIndex())
            KSpaceSolver->FromMain_AddTransducerSource();
//        //            Get_ux_sgx().AddTransducerSource(Get_u_source_index(), Get_delay_mask(), Get_transducer_source_input());
//        //
        KSpaceSolver->FromMain_Compute_duxyz();
//        //        Compute_duxyz();
//        //
//        //
        if (Parameters->Get_nonlinear_flag())
            KSpaceSolver->FromMain_Compute_rhoxyz_nonlinear();
//        //            Compute_rhoxyz_nonlinear();
        else
            KSpaceSolver->FromMain_Compute_rhoxyz_linear();
//        //            Compute_rhoxyz_linear();
//        //
//        //
//        //        // add in the source pressure term
        KSpaceSolver->FromMain_Add_p_source();
//        //        Add_p_source();
//        //
        if (Parameters->Get_nonlinear_flag()){
            KSpaceSolver->FromMain_Compute_new_p_nonlinear();
//            //            Compute_new_p_nonlinear();
        }else {
            KSpaceSolver->FromMain_Compute_new_p_linear();
//            //            Compute_new_p_linear();
        }
//        //
//        //
//        //        //-- calculate initial pressure
//        //        if ((t_index == 0) && (Parameters->Get_p0_source_flag() == 1)) Calculate_p0_source();
        if ((KSpaceSolver->GetTimeIndex() == 0) && (Parameters->Get_p0_source_flag() == 1))
            KSpaceSolver->FromMain_Calculate_p0_source();
//        
//        //
//        //
//        //
//        //        /* FIXME:
//        //         * ---------------------------------------------------------------------------------
//        //         * Here is where the injection of light should occur (i.e. the calls to the monte carlo simulation)
//        //         *
//        //         * Looking at 'StoreSensorData()' shows how to get access to the ultrasound data.
//        //         */
        
        /// --------------------------- Monte-Carlo -----------------------------------------------
        TRealMatrix *currentPressure = &(KSpaceSolver->FromMain_Get_Temp_1_RS3D());
        m_medium->Assign_current_pressure(currentPressure);
        da_boost->Run_seeded_MC_sim(m_medium, m_Laser_injection_coords);
        //float temp = currentPressure->GetElementFrom3D(100,100,58);
        //if (temp != 0) cout << "Temp = " << temp << endl;
//        //float *temp = currentPressure.GetRawData();
//        
//        //        //float *currentPressure = Get_Temp_1_RS3D().GetRawData();
//        //        
//        //        //-- store the initial pressure at the first time step --//
//        //
//        KSpaceSolver.FromMain_StoreSensorData();
//        //        StoreSensorData();
//        //
        KSpaceSolver->FromMain_PrintStatistics();
//        //        PrintStatisitcs();
//        //        
//        //
    }
    
}




void
AO_Sim::kWave_allocate_memory()
{
    assert(KSpaceSolver != NULL);
    try {
        KSpaceSolver->AllocateMemory();
    } catch (exception e){
        cout << "Failed!\nK-Wave panic: Not enough memory to run this simulation!\n";
        cout << e.what() << "\n";
        cerr << "K-Wave panic: Not enough memory to run this simulation!\n";
        cerr << e.what() << "\n";
        
        exit(EXIT_FAILURE);
    }
    
    cout << ".... done\n";

}


// Creates the grid for the monte-carlo simulation based upon the dimensions
// used in the k-Wave simulation.
//
void
AO_Sim::Create_MC_grid(TParameters * parameters)
{
    
    // Number of voxels in each axis.
    size_t Nx = parameters->GetFullDimensionSizes().X;
    size_t Ny = parameters->GetFullDimensionSizes().Y;
    size_t Nz = parameters->GetFullDimensionSizes().Z;
    
    // Size of each voxel along its respective axis.
    double dx = parameters->Get_dx();
    double dy = parameters->Get_dy();
    double dz = parameters->Get_dz();
    

    // Create the medium object that represents the grid for the
    // monte-carlo simulation based upon the dimensions for simulating
    // ultrasound.
    //
    m_medium = new Medium(Nx*dx,
                          Ny*dy,
                          Nz*dz);
    
    
    
    // Set the voxel sizes of the 'Medium'.
    this->m_medium->Set_dx(dx);
    this->m_medium->Set_dy(dy);
    this->m_medium->Set_dz(dz);
    
    // Set the number of voxels of the 'Medium'.
    this->m_medium->Set_Nx(Nx);
    this->m_medium->Set_Ny(Ny);
    this->m_medium->Set_Nz(Nz);
    
    
//    cout << Nx*dx << " "
//    << Ny*dy << " "
//    << Nz*dz << endl;
//    
//    cout << "Simulation time steps: " << parameters->Get_Nt() << endl;
    
}

