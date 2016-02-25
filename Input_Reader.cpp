//====================================================================================
//SOFTWARE:	Network of Eelectrically Conductive Nanocomposites (NECN)
//CODE FILE:	Input_Reader.cpp
//OBJECTIVE:	Reading the input data
//AUTHOR:		Fei Han; Angel Mora
//E-MAIL:			fei.han@kaust.edu.sa;  angel.mora@kaust.edu.sa
//====================================================================================

#include "Input_Reader.h"

//---------------------------------------------------------------------------
//Read data
int Input::Read_Infile(ifstream &infile)
{

	cout << "Reading input file..." << endl;
	hout << "Reading input file..." << endl;

	while(!infile.eof())
	{
		istringstream istr(Get_Line(infile));
		if(infile.eof()) break;
		string str_temp;
		istr >> str_temp;

		if(str_temp.empty()) continue;  //skip over the space or Enter key after every input item
		else if(str_temp=="Application_Name") { if(Read_application_name(app_name, infile)==0) return 0; }
		else if(str_temp=="Simulation_Parameters")	{ if(Read_simulation_parameters(simu_para, infile)==0) return 0; }
		else if(str_temp=="RVE_Geometry")	{ if(Read_rve_geometry(geom_rve, infile)==0) return 0; }
		else if(str_temp=="Nanotube_Geometry")	{ if(Read_nanotube_geo_parameters(nanotube_geo, infile)==0) return 0; }
		else if(str_temp=="Electrical_Parameters")	{ if(Read_electrical_paramters(electric_para, infile)==0) return 0; }
		else 
		{ 
			cout << "Error: the keywords \"" << str_temp << "\" is not defined!" << endl; 
			hout << "Error: the keywords \"" << str_temp << "\" is not defined!" << endl; 
			return 0; 
		}

		//the real area of cnts in the RVE
		if(geom_rve.mark&&nanotube_geo.mark)
		{
			//Calculate the real area of nanotubes
			nanotube_geo.real_area = nanotube_geo.area_fraction * geom_rve.area;
			//Define the extented RVE with one length of nanotube
			geom_rve.ex_origin.x = geom_rve.origin.x - nanotube_geo.len_max;
			geom_rve.ex_origin.y = geom_rve.origin.y - nanotube_geo.len_max;
			geom_rve.ex_len = geom_rve.len_x + 2*nanotube_geo.len_max;
			geom_rve.ey_wid = geom_rve.wid_y+ 2*nanotube_geo.len_max;
		}
	}

	cout << "Reading the keywords is finished!" << endl;
	hout << "Reading the keywords is finished!" << endl;

	if(!app_name.mark) { cout << "Attention: \"Application_Name\" will use default parameters!" << endl; hout << "Attention: \"Application_Name\" will use default parameters!" << endl; }
	if(!simu_para.mark) { cout << "Attention: \"Simulation_Parameters\" will use default parameters!" << endl; hout << "Attention: \"Simulation_Parameters\" will use default parameters!" << endl; }
	if(!geom_rve.mark) { cout << "Attention: \"RVE_Geometry\" will use default parameters!" << endl; hout << "Attention: \"RVE_Geometry\" will use default parameters!" << endl; }
	if(!nanotube_geo.mark) {	cout << "Attention: \"Nanotube_Geometry\" will use default parameters!" << endl; hout << "Attention: \"Nanotube_Geometry\" will use default parameters!" << endl; }	
	if(!electric_para.mark) {	cout << "Attention: \"Electrical_Parameters\" will use default parameters!" << endl; hout << "Attention: \"Electrical_Parameters\" will use default parameters!" << endl; }

	return 1;
}
//---------------------------------------------------------------------------
//Initialize data
int Input::Data_Initialization()
{
	//Initialize name of simulation
	app_name.keywords = "Application_Name";
	app_name.mark = false;
	app_name.str = "App_Electrical_Network_2D";

	//Initialize paramters of simulation
	simu_para.keywords = "Simulation_Parameters";
	simu_para.mark = false;
	simu_para.simu_name = "Test";
	simu_para.sample_num = 1;
	simu_para.create_read_network = "Create_Network";

	//Initialize the geometric parameters of the RVE
	geom_rve.keywords = "RVE_Geometry";
	geom_rve.mark = false;
	geom_rve.origin.x = 0.0;
	geom_rve.origin.y = 0.0;
	geom_rve.origin.flag = 0;
	geom_rve.len_x = 1.0;
	geom_rve.wid_y = 1.0;
	geom_rve.ex_origin.x = 0.0;
	geom_rve.ex_origin.y = 0.0;
	geom_rve.ex_origin.flag = 0;
	geom_rve.ex_len = 1.0;
	geom_rve.ey_wid = 1.0;
	geom_rve.area = geom_rve.len_x*geom_rve.wid_y;
	geom_rve.density = 1.0;
	geom_rve.gs_minx = 1.0;
	geom_rve.gs_miny = 1.0;
	geom_rve.win_max_x = 1.0;
	geom_rve.win_max_y = 1.0;
	geom_rve.win_min_x = 1.0;
	geom_rve.win_min_y = 1.0;
	geom_rve.win_delt_x = 1.0;
	geom_rve.win_delt_y = 1.0;

	//Initialize the geometric paramters of nanotubes
	nanotube_geo.keywords = "Nanotube_Geometry";
	nanotube_geo.mark = false;
	nanotube_geo.dir_distrib_type = "random";
	nanotube_geo.ini_pha = 0.0;
	nanotube_geo.ini_sita = 0.0;
	nanotube_geo.angle_max = 1.5707963267948966;
	nanotube_geo.step_length = 0.01;
	nanotube_geo.len_distrib_type = "uniform";
	nanotube_geo.len_min = 1.0;
	nanotube_geo.len_max = 1.0;
	nanotube_geo.rad_distrib_type = "uniform";
	nanotube_geo.rad_min = 0.005;
	nanotube_geo.rad_max = 0.005;
	nanotube_geo.criterion = "vol";
	nanotube_geo.area_fraction = 0.0;
	nanotube_geo.accum_mode = 0;
	nanotube_geo.real_area = 0.0;
	nanotube_geo.weight_fraction = 0.0;
	nanotube_geo.real_weight = 0.0;
	nanotube_geo.matrix_density = 1.06;
	nanotube_geo.linear_density = 5.8E-5;


	//Initialize electrical parameters
	electric_para.keywords = "Electrical_Parameters";
	electric_para.mark = false;
	electric_para.applied_voltage = 1.0;
	electric_para.resistivity_Ag = 0.001;

	cout << "^_^ Data initialization achieves" <<endl<<endl;
	hout << "^_^ Data initialization achieves" <<endl<<endl;

	return 1;
}
//---------------------------------------------------------------------------
//Reading the name of application case
int Input::Read_application_name(struct App_name &app_name, ifstream &infile)
{
	if(app_name.mark)
	{
		cout << "Attention: \"" << app_name.keywords << "\" has been input!" << endl;
		hout << "Attention: \"" << app_name.keywords << "\" has been input!" << endl;
		return 0;
	}
	else app_name.mark = true;

	istringstream istr(Get_Line(infile));
	istr >> app_name.str;

	return 1;
}
//---------------------------------------------------------------------------
//Reading the parameters of simulation
int Input::Read_simulation_parameters(struct Simu_para &simu_para, ifstream &infile)
{
	if(simu_para.mark)
	{
		cout << "Attention: \"" << simu_para.keywords << "\" has been input!" << endl;
		hout << "Attention: \"" << simu_para.keywords << "\" has been input!" << endl;
		return 0;
	}
	else simu_para.mark = true;

	istringstream istr0(Get_Line(infile));
	istr0 >> simu_para.simu_name;			//Read the name of simulation

	istringstream istr1(Get_Line(infile));
	istr1 >> simu_para.sample_num;			//Read the number of samples
	if(simu_para.sample_num<1)	 {	hout << "Error: the number of samples less than 1." << endl; return 0; }

	istringstream istr2(Get_Line(infile));		
	istr2 >> simu_para.create_read_network;		//Read a signal to show if create a new network or read a previouse network from a file
	if(simu_para.create_read_network!="Create_Network"&&simu_para.create_read_network!="Read_Network")
	{ hout << "Error: the 'create_read_network' is neither 'Create_Network' nor 'Read_Network'." << endl; return 0; }

	return 1;
}
//---------------------------------------------------------------------------
//Reading geometric information of the RVE
int Input::Read_rve_geometry(struct Geom_RVE &geom_rve, ifstream &infile)
{
	if(geom_rve.mark)
	{
		cout << "Attention: \"" << geom_rve.keywords << "\" has been input!" << endl;
		hout << "Attention: \"" << geom_rve.keywords << "\" has been input!" << endl;
		return 0;
	}
	else geom_rve.mark = true;

	//-----------------------------------------------------------------------------------------------------------------------------------------
	//Define the domain of RVE: the lower-left corner point of RVE and the length, width and height of RVE
	istringstream istr0(Get_Line(infile));
	istr0 >> geom_rve.origin.x >> geom_rve.origin.y;
	istr0 >> geom_rve.len_x >> geom_rve.wid_y; 
	if(geom_rve.len_x<=0||geom_rve.wid_y<=0) 
	{
		cout << "Error: the sizes of RVE should be positive!" << endl;
		hout << "Error: the sizes of RVE should be positive!" << endl;
		return 0;
	}
	geom_rve.area = geom_rve.len_x*geom_rve.wid_y; 

	//-----------------------------------------------------------------------------------------------------------------------------------------
	//Define the size range of the observation window and descrement by every step in x, y and z directions
	istringstream istr1(Get_Line(infile));
	istr1 >> geom_rve.win_max_x >> geom_rve.win_max_y; 
	istringstream istr2(Get_Line(infile));
	istr2 >> geom_rve.win_delt_x >> geom_rve.win_delt_y; 
	istringstream istr3(Get_Line(infile));
	istr3 >> geom_rve.win_min_x >> geom_rve.win_min_y; 

	if(geom_rve.win_max_x<=0.0||geom_rve.win_max_y<=0.0||
	   geom_rve.win_max_x>geom_rve.len_x||geom_rve.win_max_y>geom_rve.wid_y) 
	{
		cout << "Error: the win_max in each direction of RVE should be positive and must be smaller than the size of RVE." << endl;
		hout << "Error: the win_max in each direction of RVE should be positive and must be smaller than the size of RVE." << endl;
		return 0;
	}
	if(geom_rve.win_min_x<=0.0||geom_rve.win_min_y<=0.0||
		geom_rve.win_min_x>geom_rve.win_max_x||geom_rve.win_min_y>geom_rve.win_max_y) 
	{
		cout << "Error: the win_min in each direction of RVE should be positive and must be smaller than max." << endl;
		hout << "Error: the win_min in each direction of RVE should be positive and must be smaller than max." << endl;
		return 0;
	}
	if(geom_rve.win_delt_x<=0.0||geom_rve.win_delt_y<=0.0) 
	{
		cout << "Error: the win_delt in each direction of RVE should be positive." << endl;
		hout << "Error: the win_delt in each direction of RVE should be positive." << endl;
		return 0;
	}

	//Details: +Zero for reducing the error of division
	int num[3] = { (int)((geom_rve.win_max_x - geom_rve.win_min_x + Zero) / geom_rve.win_delt_x),
		(int)((geom_rve.win_max_y - geom_rve.win_min_y + Zero) / geom_rve.win_delt_y) };

	if(num[0]!=num[1]||num[0]!=num[2])
	{
		cout << "Error: the numbers of cutoff times are different in three directions (x, y)." << endl;
		hout << "Error: the numbers of cutoff times are different in three directions (x, y)." << endl;
		return 0;
	}
	else geom_rve.cut_num = num[0];

	//-----------------------------------------------------------------------------------------------------------------------------------------
	//Define the minimum size for background grids (looking for contact points)
	istringstream istr4(Get_Line(infile));
	istr4 >> geom_rve.gs_minx >> geom_rve.gs_miny; 
	if(geom_rve.gs_minx<=0||geom_rve.gs_miny<=0)
	{
		cout << "Error: the number of segments in each direction of RVE should be positive!" << endl;
		hout << "Error: the number of segments in each direction of RVE should be positive" << endl;
		return 0;
	}
	else if((int)(geom_rve.win_max_x/geom_rve.gs_minx)>500||
			  (int)(geom_rve.win_max_y/geom_rve.gs_miny)>500)
	{
		cout << "Error: the number of divisions in one of boundary is too big (>500), which leads to the memory problem!" << endl;
		hout << "Error: the number of divisions in one of boundary is too big (>500), which leads to the memory problem!" << endl;
		return 0;	
	}

	return 1;
}
//---------------------------------------------------------------------------
//Reading the geometric parameters of nanotubes
int Input::Read_nanotube_geo_parameters(struct Nanotube_Geo &nanotube_geo, ifstream &infile)
{
	if(nanotube_geo.mark)
	{
		cout << "Attention: \"" << nanotube_geo.keywords << "\" has been input!" << endl;
		hout << "Attention: \"" << nanotube_geo.keywords << "\" has been input!" << endl;
		return 0;
	}
	else nanotube_geo.mark = true;

	//-----------------------------------------------------------------------------------------------------------------------------------------
	//Define the initial growth direction type (random or specific) in a RVE
	istringstream istr_initial_direction(Get_Line(infile));
	istr_initial_direction >> nanotube_geo.dir_distrib_type;
	if(nanotube_geo.dir_distrib_type!="random"&&nanotube_geo.dir_distrib_type!="specific"){ hout << "Error: the direction distribution type must be either random or specific." << endl;	return 0; }
	if(nanotube_geo.dir_distrib_type=="specific")
	{
		istr_initial_direction  >> nanotube_geo.ini_sita >> nanotube_geo.ini_pha;
		if(nanotube_geo.ini_sita<0||nanotube_geo.ini_sita>PI||nanotube_geo.ini_pha<0||nanotube_geo.ini_pha>=2*PI)
		{
			hout << "Error: the specified angle is not in the acceptable range of (0, 2PI)." << endl;
			return 0;
		}
	}

	//-----------------------------------------------------------------------------------------------------------------------------------------
	//Define the normal distribution range [-omega, omega] of the growth direction
	istringstream istr_angle_range(Get_Line(infile));
	istr_angle_range >> nanotube_geo.angle_max;
	if(nanotube_geo.angle_max>0.5*PI){ hout << "Error: the specified angle is not in the acceptable range of (-PI/2, PI/2)." << endl;	 return 0; }

	//-----------------------------------------------------------------------------------------------------------------------------------------
	//Define the step length (unit: micromether) of nanotube growth
	istringstream istr_step_len(Get_Line(infile));
	istr_step_len >> nanotube_geo.step_length;
	if(nanotube_geo.step_length<=0||
	   nanotube_geo.step_length>=0.25*geom_rve.len_x||
	   nanotube_geo.step_length>=0.25*geom_rve.wid_y)
	{ hout << "Error: the step length must be positive and 0.25 times lesser than the dimension of the RVE box." << endl;	return 0; }

	//-----------------------------------------------------------------------------------------------------------------------------------------
    //Define the distribution type (uniform or normal) of the length (unit: micromether) of nanotubes and the length range (min, max) of nanotubes in a RVE
    istringstream istr_cnt_len(Get_Line(infile));
    istr_cnt_len >> nanotube_geo.len_distrib_type;
    if(nanotube_geo.len_distrib_type!="uniform"&&nanotube_geo.len_distrib_type!="normal"){ hout << "Error: the distribution of the length should be either normal or uniform." << endl;	return 0; }
    istr_cnt_len >> nanotube_geo.len_min >> nanotube_geo.len_max;
    if(nanotube_geo.len_min<0||nanotube_geo.len_max<0||nanotube_geo.len_max<nanotube_geo.len_min){ hout << "Error: the length must be non-negative and min must be smaller than max." << endl; return 0; }

	//-----------------------------------------------------------------------------------------------------------------------------------------
	//Define the distribution type (uniform or normal) of the radius (unit: micromether) of nanotubes and the radius range (min, max) of nanotubes in a RVE
    istringstream istr_cnt_rad(Get_Line(infile));
    istr_cnt_rad >> nanotube_geo.rad_distrib_type;
    if(nanotube_geo.rad_distrib_type!="uniform"&&nanotube_geo.rad_distrib_type!="normal"){ hout << "Error: the distribution of the radius should be either normal or uniform." << endl;	return 0; }
    istr_cnt_rad >> nanotube_geo.rad_min >> nanotube_geo.rad_max;
    if(nanotube_geo.rad_min<0||nanotube_geo.rad_max<0||nanotube_geo.rad_max<nanotube_geo.rad_min||
	   nanotube_geo.rad_min>3*nanotube_geo.step_length||nanotube_geo.rad_max>0.05*nanotube_geo.len_min)
	{ hout << "Error: the radius must be non-negative, min must be smaller than max, min must be smaller than 3*step_length and max must be smaller than 0.05*len_min." << endl; return 0; }

	//-----------------------------------------------------------------------------------------------------------------------------------------
	//Define the volume or weight fraction of nanotubes in the RVE
	istringstream istr_cnt_vol(Get_Line(infile));
	istr_cnt_vol >> nanotube_geo.criterion;
	if(nanotube_geo.criterion=="area")
	{
		istr_cnt_vol >> nanotube_geo.area_fraction;
		if(nanotube_geo.area_fraction>1||nanotube_geo.area_fraction<0){ hout << "Error: the area fraction must be between 0 and 1." << endl; return 0; }
		hout << "    The area fraction is "<< nanotube_geo.area_fraction << endl;
        
		istr_cnt_vol >> nanotube_geo.accum_mode;
		if(nanotube_geo.accum_mode<0&&nanotube_geo.accum_mode>2){ hout <<"Error: the mode of accumulation should be between 0 and 2." << endl; return 0; }
        
		//The total volume of the nanotube network
		nanotube_geo.real_area = nanotube_geo.area_fraction*geom_rve.area;
	}
	else if(nanotube_geo.criterion=="wt")
	{
		istr_cnt_vol >> nanotube_geo.weight_fraction;
		if(nanotube_geo.weight_fraction>1||nanotube_geo.weight_fraction<0){ hout << "Error: the area fraction must be between 0 and 1." << endl; return 0; }
		hout << "    The weight fraction is " << nanotube_geo.weight_fraction << endl;
        
		istr_cnt_vol >> nanotube_geo.accum_mode;
		if(nanotube_geo.accum_mode<0&&nanotube_geo.accum_mode>2){ hout <<"Error: the mode of accumulation should be between 0 and 2." << endl; return 0;  }
        
		istr_cnt_vol >> nanotube_geo.linear_density;		//Read the linear density of a nanotube
		if(nanotube_geo.linear_density<0){ hout << "Error: the linear density of a nanotube should be non-nagetive." << endl; return 0; }
		
		istr_cnt_vol >> geom_rve.density;	 //Read the density of RVE. Here we ignore the volume of nantubes, so the density of RVE actually approximates to the density of matix
		if(geom_rve.density<0){ hout << "Error: the density of RVE should be non-nagetive." << endl; return 0; }
		if(nanotube_geo.linear_density>=nanotube_geo.matrix_density){ hout << "Error: the density of matrix or the linear density of a nanotube is wrong." << endl; return 0; }
        
		//The real weight of nanotubes
		nanotube_geo.real_weight = nanotube_geo.weight_fraction*geom_rve.area*geom_rve.density;
	}
	else { hout << "Error: the criterian of generation is neither 'area' nor 'wt'." << endl; return 0; }

	return 1;
}

//---------------------------------------------------------------------------
//Reading current information
int Input::Read_electrical_paramters(struct Electric_para &electric_para, ifstream &infile)
{
	if(electric_para.mark)
	{
		cout << "Attention: \"" << electric_para.keywords << "\" has been input!" << endl;
		hout << "Attention: \"" << electric_para.keywords << "\" has been input!" << endl;
		return 0;
	}
	else electric_para.mark = true;

	istringstream istr0(Get_Line(infile));
	istr0 >> electric_para.applied_voltage;		

	istringstream istr1(Get_Line(infile));
	istr1 >> electric_para.resistivity_Ag;

	return 1;
}
//---------------------------------------------------------------------------
//Read the input data in a whole line (to skip over the comment line starting with a '%')
string Input::Get_Line(ifstream &infile)const
{
	string s;
	//Read the input data in a whole line
	getline(infile,s);
	//to skip over the comment line starting with a '%'
	while(!infile.eof() && s.substr(0,1)=="%")
		getline(infile,s);
	return s;
}
//===========================================================================
