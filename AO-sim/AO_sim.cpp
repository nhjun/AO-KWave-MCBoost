//
//  AO_sim.cpp
//  K_Wave_C
//
//  Created by jacob on 12/3/12.
//  Copyright (c) 2012 BMPI. All rights reserved.
//

#include <exception>

#include "AO_sim.h"
#include <MC-Boost/logger.h>





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

	/// Set default value.
	MC_time_step = 0.0f;

}



AO_Sim::~AO_Sim()
{
    cout << "\n\nAO_Sim:: destructor\n";
    if (da_boost)
        delete da_boost;

    if (KSpaceSolver)
        delete KSpaceSolver;

//    if (m_medium)
//        delete m_medium;

}



/// Run the actual acousto-optic simulation.
void
AO_Sim::Run_acousto_optics_sim(TParameters * Parameters,
							   bool sim_displacement,
							   bool sim_refractive_grad)
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
    KSpaceSolver->FromAO_sim_PrintOutputHeader();
    KSpaceSolver->IterationTimeStart();
	size_t k_wave_Nt = Parameters->Get_Nt();
	k_wave_Nt = 901;
    for (KSpaceSolver->SetTimeIndex(0); KSpaceSolver->GetTimeIndex() < k_wave_Nt; KSpaceSolver->IncrementTimeIndex()){

        cout << ".......... Running k-Wave ........... ("
             << KSpaceSolver->GetTimeIndex() << " of "
             << k_wave_Nt << '\n';
            //<< Parameters->Get_Nt() << ")\n";

        KSpaceSolver->FromAO_sim_Compute_uxyz();

	// add in the velocity u source term
        KSpaceSolver->FromAO_sim_Add_u_source();

	// add in the transducer source term (t = t1) to ux
        if (Parameters->Get_transducer_source_flag() > KSpaceSolver->GetTimeIndex())
            KSpaceSolver->FromAO_sim_AddTransducerSource();


        KSpaceSolver->FromAO_sim_Compute_duxyz();


        if (Parameters->Get_nonlinear_flag())
            KSpaceSolver->FromAO_sim_Compute_rhoxyz_nonlinear();
        else
            KSpaceSolver->FromAO_sim_Compute_rhoxyz_linear();


	// add in the source pressure term
        KSpaceSolver->FromAO_sim_Add_p_source();

        if (Parameters->Get_nonlinear_flag()){
            KSpaceSolver->FromAO_sim_Compute_new_p_nonlinear();
        }else {
            KSpaceSolver->FromAO_sim_Compute_new_p_linear();

        }


	//-- calculate initial pressure
        if ((KSpaceSolver->GetTimeIndex() == 0) && (Parameters->Get_p0_source_flag() == 1))
            KSpaceSolver->FromAO_sim_Calculate_p0_source();



        //-- store the initial pressure at the first time step --//
        /// Here is where the computation for the refractive index, displacements,
        /// and all the values kWave needs, takes place if commandline flags for saving data hvave been set.
        KSpaceSolver->FromAO_sim_StoreSensorData();


        KSpaceSolver->FromAO_sim_PrintStatistics();



        /// --------------------------- Begin Monte-Carlo Simulation ------------------------------------------------------
        ///
		/// If the flag for simulating the influence of refractive index changes or displacements is set
		/// we need to grab the appropriate data from k-Wave and build grids for
		/// photon propagation.
		if (sim_refractive_grad || sim_displacement )
		{
            TRealMatrix *refractive_x = NULL;
            TRealMatrix *refractive_y = NULL;
            TRealMatrix *refractive_z = NULL;

            TRealMatrix *disp_x = NULL;
            TRealMatrix *disp_y = NULL;
            TRealMatrix *disp_z = NULL;

            if (sim_refractive_grad)
            {
                /// Check if we are saving the refracive index data in KSpaceSolver, which is based on a commandline flag.
                /// If we are, we don't want to run the computation from here, as it will already
                /// be done in KSpaceSolver if we are saving the data.  This eliminates the possiblity
                /// of computing something twice, which could happen when monte-carlo needs the data
                /// and the flag has been set to save data in KSpaceSolver.
                if (!(Parameters->IsStore_refractive_x() ||
                      Parameters->IsStore_refractive_y() ||
                      Parameters->IsStore_refractive_z()))
                {
                    /// This is executed when the AO simulation needs refractive index changes, but
                    /// the commandline flags were not made to save the refractive index data, as
                    /// discussed above.
                    KSpaceSolver->FromAO_sim_compute_refractive_index();
                }

                refractive_x = &(KSpaceSolver->FromAO_sim_Get_refractive_x());
                refractive_y = &(KSpaceSolver->FromAO_sim_Get_refractive_y());
                refractive_z = &(KSpaceSolver->FromAO_sim_Get_refractive_z());
            }

            if (sim_displacement)
            {
                /// Same reasoning as above.
                if (!(Parameters->IsStore_disp_x() ||
                      Parameters->IsStore_disp_y() ||
                      Parameters->IsStore_disp_z()))
                {
                    KSpaceSolver->FromAO_sim_compute_displacement();
                }

                disp_x = &(KSpaceSolver->FromAO_sim_Get_disp_x());
                disp_y = &(KSpaceSolver->FromAO_sim_Get_disp_y());
                disp_z = &(KSpaceSolver->FromAO_sim_Get_disp_z());
            }

			/// NOTE:
			/// - The refractive index and displacement values need to be updated every time step
            ///   within kWave, otherwise what is assigned here is invalid.
			if (sim_refractive_grad &&
				sim_displacement)
			{
				/// Create a refractive map based upon the pressure and density changes at this time step.
        		m_medium->Create_refractive_map(refractive_x,
                                                refractive_y,
                                                refractive_z);

				/// Create a displacment map based upon the pressure at this time step.
        		m_medium->Create_displacement_map(disp_x,
                                                  disp_y,
                                                  disp_z);
			}
			else if (sim_refractive_grad)
			{
                /// Create a refractive map based upon the pressure and density changes at this time step.
        		m_medium->Create_refractive_map(refractive_x,
                                                refractive_y,
                                                refractive_z);
			}
			/// Similar to above (i.e. sim_refractive_grad)
			else if (sim_displacement)
			{
                /// Create a displacment map based upon the pressure at this time step.
        		m_medium->Create_displacement_map(disp_x,
                                                  disp_y,
                                                  disp_z);
        	}


			/// Decide what to simulate (refractive gradient, displacement) based on whether those objects exist.
			m_medium->kwave.nmap == NULL ? da_boost->Simulate_refractive_gradient(false) :
                                           da_boost->Simulate_refractive_gradient(true);

			/// Decide what to simulate (refractive gradient, displacement) based on whether those objects exist.
        	m_medium->kwave.dmap == NULL ? da_boost->Simulate_displacement(false) :
                                           da_boost->Simulate_displacement(true);

        	/// Not saving seeds, so set to false.
        	da_boost->Save_RNG_Seeds(false);


        	/// Only run the MC-sim after ultrasound has propagated a certain distance (or time).
			/// Similar to the stroboscopic experiments.
        	static size_t cnt = MC_time_step/Parameters->Get_dt();

        	if ((((KSpaceSolver->GetTimeIndex() % cnt) == 0) && (KSpaceSolver->GetTimeIndex() >= 700)
                 && (KSpaceSolver->GetTimeIndex() <= 900)) || KSpaceSolver->GetTimeIndex() == 1)
        	{
            	cout << ".......... Running MC-Boost ......... ";
            	cout << "(time: " << KSpaceSolver->GetTimeIndex()*Parameters->Get_dt() << ")\n";
            	da_boost->Run_seeded_MC_sim_timestep(m_medium,
            	                                     m_Laser_injection_coords,
            	                                     KSpaceSolver->GetTimeIndex());
        	}
		} // END if (sim_refractive_grad || sim_displacement )

        ///
        /// --------------------------- End Monte-Carlo Simulation ------------------------------------------------------

    }


    KSpaceSolver->IterationTimeStop();

    cout << "-------------------------------------------------------------\n";
    cout << "Elapsed time: " << KSpaceSolver->GetSimulationTime() << "\n";
    cout << "-------------------------------------------------------------\n";

    cout << "Post-processing phase......."; cout.flush();
    KSpaceSolver->PostProcessingTimeStart();


    KSpaceSolver->FromAO_sim_PostProcessing();
    KSpaceSolver->PostProcessingTimeStop();
    cout << "Done \n";
    cout << "Elapsed time: " << KSpaceSolver->GetPostProcessingTime() << "\n";


    KSpaceSolver->FromAO_sim_WriteOutputDataInfo();

    Parameters->HDF5_OutputFile.Close();


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
    // Number of voxels in each axis, with the PML taken into account.
    // NOTE:
    // - It is assumed that the PML is always included INSIDE the computational
    //   medium (i.e. it is not expanded), as defined in the matlab script, and that there is a small
    //   cushion (5 voxels along US transmission axis) to ensure the probe sits away from the PML.
    //   These things are taken into account below in the monte-carlo simulation grid.
    // - As an example, if the x-axis is 512 elements and the X_PML is 20 elements,
    //   the US probe should sit at X_PML+OFFSET in the k-Wave computation medium, which
    //   means the monte-carlo medium should reflect these changes because US data in the PML
    //   is heavily distorted, as it should be, and we don't want photons moving through those
    //   regions.
    size_t OFFSET = 5; /// An offset to apply to the transmission axis PML (i.e. x_pml).
    size_t x_pml_offset = parameters->Get_pml_x_size()+OFFSET;
    size_t y_pml_offset = parameters->Get_pml_y_size();
    size_t z_pml_offset = parameters->Get_pml_z_size();

    size_t Nx = parameters->GetFullDimensionSizes().X - 2*x_pml_offset;
    size_t Ny = parameters->GetFullDimensionSizes().Y - 2*y_pml_offset;
    size_t Nz = parameters->GetFullDimensionSizes().Z - 2*z_pml_offset;

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

    /// Set the size of the PML's used in the k-Wave simulation.
    this->m_medium->Set_X_PML(x_pml_offset);
    this->m_medium->Set_Y_PML(y_pml_offset);
    this->m_medium->Set_Z_PML(z_pml_offset);


    /// Set the size of the sensor mask used in the kWave simulation.  This is the number
    /// of elements in the TRealMatrix that recorded data.
    // const size_t  sensor_size = Get_sensor_mask_ind().GetTotalElementCount();
    this->m_medium->kwave.sensor_mask_index_size    = parameters->Get_sensor_mask_index_size();

    /// Set the time-step used in k-Wave.
    this->m_medium->kwave.dt = parameters->Get_dt();
    //this->m_medium->kwave.speed_of_sound            = parameters->Get_c0_scalar();


}




