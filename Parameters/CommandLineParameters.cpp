/**
 * @file        CommandLineParameters.cpp
 * @author      Jiri Jaros              \n
 *              CECS, ANU, Australia    \n
 *              jiri.jaros@anu.edu.au
 * @brief       The implementation file containing the command line parameters
 *
 * @version     kspaceFirstOrder3D 2.13
 * @date        29 August 2012, 11:25 (created) \n
 *              11 October 2012, 17:05 (revised)
 *
 *
 *
 * @section License
 * This file is part of the C++ extension of the k-Wave Toolbox (http://www.k-wave.org).\n
 * Copyright (C) 2012 Jiri Jaros and Bradley Treeby
 *
 * This file is part of k-Wave. k-Wave is free software: you can redistribute it
 * and/or modify it under the terms of the GNU Lesser General Public License as
 * published by the Free Software Foundation, either version 3 of the License,
 * or (at your option) any later version.
 *
 * k-Wave is distributed in the hope that it will be useful, but
 * WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
 * See the GNU Lesser General Public License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public License
 * along with k-Wave. If not, see <http://www.gnu.org/licenses/>.
 */

#include <stdio.h>
#include <getopt.h>
#include <string.h>
#include <omp.h>

#include <Parameters/CommandLineParameters.h>

#include <Utils/ErrorMessages.h>
//----------------------------------------------------------------------------//
//---------------------------- Constants -------------------------------------//
//----------------------------------------------------------------------------//

//----------------------------------------------------------------------------//
//----------------------------- Public   -------------------------------------//
//----------------------------------------------------------------------------//



/**
 * Constructor
 */
TCommandLineParameters::TCommandLineParameters() :
        InputFileName(""), OutputFileName (""),
        NumberOfThreads(omp_get_num_procs()),
        VerboseInterval(DefaultVerboseInterval), CompressionLevel (DefaultCompressionLevel),
        BenchmarkFlag (false), BenchmarkTimeStepsCount(0),
        PrintVersion (false),
        Store_p_raw(false), Store_p_rms(false), Store_p_max(false), Store_p_final(false),
        Store_u_raw(false), Store_u_rms(false), Store_u_max(false), Store_u_final(false),
        Store_I_avg(false), Store_I_max(false),
        /// ------------------- JWJS ----------------------
        Store_refractive_total(false), Store_refractive_x(false), Store_refractive_y(false), Store_refractive_z(false),
        Store_disp_x(false), Store_disp_y(false), Store_disp_z(false),
        /// -------------------------
        StartTimeStep(0), EndTimeStep(-1)
{

}// end of TCommandLineParameters
//------------------------------------------------------------------------------

/**
 * Print usage and exit.
 */
void TCommandLineParameters::PrintUsageAndExit(){


 printf("---------------------------------- Usage ---------------------------------\n");
 printf("Mandatory parameters:\n");
 printf("  -i <input_file_name>            : HDF5 input file\n");
 printf("  -o <output_file_name>           : HDF5 output file\n");
 printf("\n");
 printf("Optional parameters: \n");
 printf("  -t <num_threads>                : Number of CPU threads\n");
 printf("                                      (default = %d)\n",omp_get_num_procs());
 printf("  -r <interval_in_%%>              : Progress print interval\n");
 printf("                                      (default = %d%%)\n",DefaultVerboseInterval);
 printf("  -c <comp_level>                 : Output file compression level <0,9>\n");
 printf("                                      (default = %d)\n",DefaultCompressionLevel );
 printf("  --benchmark <steps>             : Run a specified number of time steps\n");
 printf("\n");
 printf("  -h                              : Print help\n");
 printf("  --help                          : Print help\n");
 printf("  --version                       : Print version\n");
 printf("\n");
 printf("Output flags:\n");
 printf("  -p                              : Store acoustic pressure \n");
 printf("                                      (default if nothing else is on)\n");
 printf("                                      (the same as --p_raw)\n");
 printf("  --p_raw                         : Store raw time series of p (default)\n");
 printf("  --p_rms                         : Store rms of p\n");
 printf("  --p_max                         : Store max of p\n");
 printf("  --p_final                       : Store final pressure field \n");
 printf("\n");
 printf("  -u                              : Store ux, uy, uz\n");
 printf("                                      (the same as --u_raw)\n");
 printf("  --u_raw                         : Store raw time series of ux, uy, uz\n");
 printf("  --u_rms                         : Store rms of ux, uy, uz\n");
 printf("  --u_max                         : Store max of ux, uy, uz\n");
 printf("  --u_final                       : Store final acoustic velocity\n");
 printf("\n");
 printf("  -I                              : Store intensity\n");
 printf("                                      (the same as --I_avg)\n");
 printf("  --I_avg                         : Store avg of intensity\n");
 printf("  --I_max                         : Store max of intensity\n");
 printf("\n");
 /// --------------------- JWJS ---------------------------------------------------------------------
 printf("  -n                              : Store index of refraction\n");
 printf("                                       (all axial components nx, ny, nz)\n");
 printf("  --refractive_total              : Store the norm of the index of refraction\n");
 printf("  --refractive_x                  : Store the x-component of the index of refraction\n");
 printf("  --refractive_y                  : Store the y-component of the index of refraction\n");
 printf("  --refractive_z                  : Store the z-component of the index of refraction\n");
 printf("\n");
 printf("  -d                              : Store displacements\n");
 printf("                                       (all axial components disp_x, disp_y, disp_z\n");
 printf("  --disp_x                        : Store displacement along x-axis\n");
 printf("  --disp_y                        : Store displacement along y-axis\n");
 printf("  --disp_z                        : Store displacement along z-axis\n");
 printf("\n");
 printf("  -e <timestep>                    : Time step when data collection ends\n");
 /// --------------------------------------------------------
 printf("\n");
 printf("  -s <timestep>                   : Time step when data collection begins\n");
 printf("                                      (default = 1)\n");
 printf("--------------------------------------------------------------------------\n");
 printf("\n");


 exit(EXIT_FAILURE);

}// end of PrintUsageAndExit
//------------------------------------------------------------------------------

/**
 * Print setup.
 */
void TCommandLineParameters::PrintSetup(){

    printf("List of enabled parameters:\n");

    printf("  Input  file           %s\n",InputFileName.c_str());
    printf("  Output file           %s\n",OutputFileName.c_str());
    printf("\n");
    printf("  Number of threads     %d\n", NumberOfThreads);
    printf("  Verbose interval[%%]  %d\n", VerboseInterval);
    printf("  Compression level     %d\n", CompressionLevel);
    printf("\n");
    printf("  Benchmark flag        %d\n", BenchmarkFlag);
    printf("  Benchmark time steps  %d\n", BenchmarkTimeStepsCount);
    printf("\n");
    printf("  Store p_raw           %d\n", Store_p_raw);
    printf("  Store p_rms           %d\n", Store_p_rms);
    printf("  Store p_max           %d\n", Store_p_max);
    printf("  Store p_final         %d\n", Store_p_final);
    printf("\n");
    printf("  Store u_raw           %d\n", Store_u_raw);
    printf("  Store u_rms           %d\n", Store_u_rms);
    printf("  Store u_max           %d\n", Store_u_max);
    printf("  Store u_max           %d\n", Store_u_final);
    printf("\n");
    printf("  Store I_avg           %d\n", Store_I_avg);
    printf("  Store I_max           %d\n", Store_I_max);
    printf("\n");
    printf("  Store refractive_total         %d\n", Store_refractive_total);
    printf("  Store refractive_x              %d\n", Store_refractive_x);
    printf("  Store refractive_y              %d\n", Store_refractive_y);
    printf("  Store refractive_z              %d\n", Store_refractive_z);
    printf("\n");
    printf("  Store disp_x          %d\n", Store_disp_x);
    printf("  Store disp_y          %d\n", Store_disp_y);
    printf("  Store disp_z          %d\n", Store_disp_z);
    printf("\n");
    printf("  Collection begins at  %d\n", StartTimeStep+1);
    printf("\n");
    printf("  Collection ends at %d\n", EndTimeStep);


}// end of PrintSetup
//------------------------------------------------------------------------------

/**
 * Parse command line.
 * @param [in, out] argc
 * @param [in, out] argv
 */
void TCommandLineParameters::ParseCommandLine(int argc, char** argv){

   char c;
   int longIndex;
   const char * shortOpts = "i:o:v:c:t:s:e:puIhnd:";

   const struct option longOpts[] = {
        { "benchmark", required_argument , NULL, 0},
        { "help", no_argument, NULL, 'h' },
        { "version", no_argument, NULL, 0 },

        { "p_raw", no_argument, NULL, 'p' },
        { "p_rms", no_argument, NULL, 0 },
        { "p_max", no_argument, NULL, 0 },
        { "p_final", no_argument, NULL, 0 },

        { "u_raw", no_argument, NULL, 'u' },
        { "u_rms", no_argument, NULL, 0 },
        { "u_max", no_argument, NULL, 0 },
        { "u_final", no_argument, NULL, 0 },

        { "I_avg", no_argument, NULL, 'I' },
        { "I_max", no_argument, NULL, 0 },

        /// ---------------- JWJS -----------------------
        //{ "n", no_argument, NULL, 'n'},
        { "refractive_total", no_argument, NULL, 0},
        { "refractive_x", no_argument, NULL, 0},
        { "refractive_y", no_argument, NULL, 0},
        { "refractive_z", no_argument, NULL, 0},

        //{ "d", no_argument, NULL, 'd'},
        { "disp_x", no_argument, NULL, 0},
        { "disp_y", no_argument, NULL, 0},
        { "disp_z", no_argument, NULL, 0},

        { "s", required_argument, NULL, 's'},
        { "e", required_argument, NULL, 'e'},
        /// ----------------------

        { NULL, no_argument, NULL, 0 }
    };


   // Short parameters //
   while ((c = getopt_long (argc, argv, shortOpts, longOpts, &longIndex )) != -1){
       switch (c){

          case 'i':{
             InputFileName = optarg;
             break;
          }
          case 'o':{
             OutputFileName = optarg;
             break;
          }

          case 'v': {
              if ((optarg == NULL) || (atoi(optarg) <= 0)) {
                  fprintf(stderr,"%s", CommandlineParameters_ERR_FMT_NoVerboseIntreval);
                  PrintUsageAndExit();
              }else {
                  VerboseInterval = atoi(optarg);
              }

              break;
          }

          case 't':{
              if ((optarg == NULL) || (atoi(optarg) <= 0)) {
                  fprintf(stderr,"%s", CommandlineParameters_ERR_FMT_NoThreadNumbers);
                  PrintUsageAndExit();
              }else {
                NumberOfThreads = atoi(optarg);
              }

              break;
          }

          case 'c':{
               if ((optarg == NULL) || (atoi(optarg) < 0) || atoi(optarg) > 9) {
                  fprintf(stderr,"%s", CommandlineParameters_ERR_FMT_NoCompressionLevel);
                  PrintUsageAndExit();
               } else {
                   CompressionLevel = atoi(optarg);
               }

             break;
          }

          case 'p':{
             Store_p_raw = true;
             break;
          }

          case 'u':{
             Store_u_raw = true;
             break;
          }

          case 'I':{
             Store_I_avg = true;
             break;
          }

        /// ------------------ JWJS ------------------------------
          case 'n':{
             Store_refractive_x     = true;
             Store_refractive_y     = true;
             Store_refractive_z     = true;
             break;
          }

          case 'd':{
             Store_disp_x = true;
             Store_disp_y = true;
             Store_disp_z = true;
             break;
          }
        /// ------------------------

          case 'h':{

             PrintUsageAndExit();
             break;
          }

          case 's':{
               if ((optarg == NULL) || (atoi(optarg) < 1)) {
                  fprintf(stderr,"%s", CommandlineParameters_ERR_FMT_NoStartTimestep);
                  PrintUsageAndExit();
               }
               StartTimeStep = atoi(optarg) - 1;

             break;
          }
        /// ----------------------- JWJS ------------------------------
           case 'e':{
               if ((optarg == NULL) || (atoi(optarg) < 1)) {
                   fprintf(stderr,"%s", CommandlineParameters_ERR_FMT_NoEndTimestep);
                   PrintUsageAndExit();
               }
               EndTimeStep = atoi(optarg) - 1;

               break;

           }
           /// ------------------------

           case 0:{   /* long option without a short arg */
                if( strcmp( "benchmark", longOpts[longIndex].name ) == 0 ) {
                     BenchmarkFlag = true;
                     if ((optarg == NULL) || (atoi(optarg) <= 0)) {
                        fprintf(stderr,"%s", CommandlineParameters_ERR_FMT_NoBenchmarkTimeStepCount);
                        PrintUsageAndExit();
                      }else {
                         BenchmarkTimeStepsCount = atoi(optarg);
                     }

                }else

                if( strcmp( "version", longOpts[longIndex].name ) == 0 ) {
                    PrintVersion = true;
                    return;
                } else

                if( strcmp( "p_rms", longOpts[longIndex].name ) == 0 ) {
                    Store_p_rms = true;
                } else
                if( strcmp( "p_max", longOpts[longIndex].name ) == 0 ) {
                    Store_p_max = true;
                } else
                if( strcmp( "p_final", longOpts[longIndex].name ) == 0 ) {
                    Store_p_final = true;
                } else

                if( strcmp( "u_rms", longOpts[longIndex].name ) == 0 ) {
                    Store_u_rms = true;
                } else
                if( strcmp( "u_max", longOpts[longIndex].name ) == 0 ) {
                    Store_u_max = true;
                } else
                if( strcmp( "u_final", longOpts[longIndex].name ) == 0 ) {
                    Store_u_final = true;
                } else
                if( strcmp( "I_max", longOpts[longIndex].name ) == 0 ) {
                    Store_I_max = true;
                } else
               /// -------------------- JWJS -----------------------------
                if( strcmp( "refractive_total", longOpts[longIndex].name ) == 0) {
                    Store_refractive_total = true;
                } else
                if( strcmp( "refractive_x", longOpts[longIndex].name ) == 0) {
                    Store_refractive_x = true;
                } else
                if( strcmp( "refractive_y", longOpts[longIndex].name ) == 0) {
                    Store_refractive_y = true;
                } else
                if( strcmp( "refractive_z", longOpts[longIndex].name ) == 0) {
                    Store_refractive_z = true;
                } else
                if( strcmp( "disp_x", longOpts[longIndex].name ) == 0) {
                    Store_disp_x = true;
                } else
                if( strcmp( "disp_y", longOpts[longIndex].name ) == 0) {
                    Store_disp_y = true;
                }else
                if( strcmp( "disp_z", longOpts[longIndex].name ) == 0) {
                    Store_disp_z = true;
                }
               /// ---------------------------
                else {
                    PrintUsageAndExit();
                }
                break;
           }
          default:{
               PrintUsageAndExit();
          }
       }
   }


    //-- Post checks --//


   if (InputFileName == "") {
       fprintf(stderr,"%s",CommandlineParameters_ERR_FMT_NoInputFile);
       PrintUsageAndExit();
   }


   if (OutputFileName == "") {
       fprintf(stderr,"%s",CommandlineParameters_ERR_FMT_NoOutputFile);
       PrintUsageAndExit();
   }


   if (!(Store_p_raw || Store_p_rms || Store_p_max || Store_p_final ||
         Store_u_raw || Store_u_rms || Store_u_max || Store_u_final ||
         Store_I_avg || Store_I_max ||
       /// ------------------------- JWJS ------------------------------------
         Store_refractive_total || Store_refractive_x || Store_refractive_y || Store_refractive_z ||
         Store_disp_x || Store_disp_y || Store_disp_z))
       /// -------------------------------
   {
            Store_p_raw = true;
   }


}// end of ParseCommandLine
//------------------------------------------------------------------------------
