//====================================================================================
//SOFTWARE:	Network of Eelectrically Conductive Nanocomposites (NECN)
//CODE FILE:	Input_Reader.h
//OBJECTIVE:	Reading the input data
//AUTHOR:		Fei Han; Angel Mora
//E-MAIL:			fei.han@kaust.edu.sa	;	angel.mora@kaust.edu.sa
//====================================================================================

#ifndef INPUTREADER_H
#define INPUTREADER_H

#include<iomanip>
#include<iostream>
#include<fstream>
#include<sstream>
#include<vector>
#include<string>
#include "time.h"
#include "Hns.h"
using namespace hns;

#include "Gauss.h"
#include "Geometry_3D.h"

const double PI = 3.1415926535897932;

//---------------------------------------------------------------------------
//Name of application case
struct App_name{
			string keywords;
			bool mark;
			string str;
		};
//Name of simulation
struct Simu_para{
			string keywords;
			bool mark;
			string simu_name;
			string create_read_network;
			int sample_num;
		};
//The geometry of the RVE
struct Geom_RVE{
			string keywords;
			bool mark;
			Point_3D origin;								//Define an origin point for a RVE
			double len_x, wid_y, hei_z;				//Define length, width and height for an  RVE for generation with an acurrate control 
			Point_3D ex_origin;							//Define an origin point for an extended RVE to generate network with an acurrate control  
			double ex_len, ey_wid, ez_hei;			//Define length, width for an extended RVE for generation with an acurrate control  
			double volume, crosec_area;				//Define the volume and the cross-sectional area
			double density;
			double gs_minx, gs_miny;					//Define the minimum size for background grids (looking for contact points)
			double win_max_x, win_max_y;		//Define the size range of the cutoff window and descrement by every step in x, y directions
			double win_min_x, win_min_y;
			double win_delt_x, win_delt_y;
			int cut_num;										//Define the number of cutoff times (0: the maxmum size, n: the maxmum size - n*step_length(delta), n>=1)
		};
//The nanowire parameters in a network
struct Nanowire_Geo{
			string keywords;
			bool mark;
			string criterion;						//Define the area or weight fraction of nanowires in the RVE: vol, wt, nwt
			string dir_distrib_type;			//Define the initial growth direction type (random or specific) in a RVE
			string len_distrib_type;			//Define the distribution type (uniform or normal) of the length (unit: micromether) of nanowires
			string rad_distrib_type;			//Define the distribution type (uniform or normal) of the radius (unit: micromether) of nanowires
			double step_length;				//Define the step length of nanotube growth
			double angle_min, angle_max;	//Define the angle range of nanowires which is the bound for a given orientational distribution type
			double len_min, len_max;		//Define the length range (min, max) of nanowires
			double rad_min, rad_max;		//Define the radius range (min,max) of nanowires
			double area_fraction;				//Define the area fraction of nanowires
			int accum_mode;					//Define the mode of accumulator (0: no accumu; 1: linear accum as sample number; 2: square exponential accum as sample (number-1))
			double real_area;					// Define the real area of nanowires
			double weight_fraction;			//Define the weight fraction of nanowires
			double real_weight;				//Define the real weight of nanowires
			double linear_density;			//Define the linear density of nanowires
			double matrix_density;			//Define the density of matrix
		};
//The cutoff distances
struct Cutoff_dist{
			string keywords;
			bool mark;
			double van_der_Waals_dist;
			double tunneling_dist;
		};
//The electrical parameters
struct Electric_para{
			string keywords;
			bool mark;
			double applied_voltage;			//Define the magnitude of the applied voltage
			double resistivity_NW;			//Define the resistivity value of the nanowire
		};

//---------------------------------------------------------------------------
class Input 
{
	public:
		//Data members
		struct App_name app_name;
		struct Simu_para simu_para;
		struct Geom_RVE geom_rve;
		struct Nanowire_Geo nanowire_geo;
		struct Cutoff_dist cutoff_dist;
		struct Electric_para electric_para;

		//Constructor
		Input(){};  

		//Member functions
		int Data_Initialization();								//Initialize data
		int Read_Infile(ifstream &infile);				//Read data
		string Get_Line(ifstream &infile)const;		//Read the input data in a whole line (to skip over the comment line starting with a '%')
private:
		//Member functions
		int Read_application_name(struct App_name &app_name, ifstream &infile);
		int Read_simulation_parameters(struct Simu_para &simu_para, ifstream &infile);
		int Read_rve_geometry(struct Geom_RVE &geom_rve, ifstream &infile);
		int Read_nanowire_geo_parameters(struct Nanowire_Geo &nanowire_geo, ifstream &infile);
		int Read_cutoff_distances(struct Cutoff_dist &cutoff_dist, ifstream &infile);
		int Read_electrical_paramters(struct Electric_para &electric_para, ifstream &infile);
};
//---------------------------------------------------------------------------
#endif
//===========================================================================
