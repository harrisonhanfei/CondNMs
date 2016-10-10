//====================================================================================
//SOFTWARE:	Network of Eelectrically Conductive Nanocomposites (NECN)
//CODE FILE:	Tecplot_Export.h
//OBJECTIVE:	To export the 3D geometric images through Tecplot data files
//AUTHOR:		Fei Han; Angel Mora
//E-MAIL:			fei.han@kaust.edu.sa	;	angel.mora@kaust.edu.sa
//====================================================================================

#include "Tecplot_Export.h"

//---------------------------------------------------------------------------
//The geometric structure of CNT network (by threads in Tecplot) in a cuboid
int Tecplot_Export::Export_network_threads(const struct cuboid &cub, const vector<vector<Point_3D> > &cnts_points)const
{
	ofstream otec("CNT_Wires.dat");
	otec << "TITLE = CNT_Wires" << endl;
	otec << "VARIABLES = X, Y, Z" << endl;

	//---------------------------------------------------------------------------
	//Export a 3D cuboid
	if(Export_cuboid(otec, cub)==0) return 0;
	
	//---------------------------------------------------------------------------
	//Export 3D nanotube threads
	if(Export_nano_threads(otec, cub, cnts_points)==0) return 0;

	//---------------------------------------------------------------------------
	otec.close();

	return 1;
}
//---------------------------------------------------------------------------
//The geometric structure of CNT network (by threads in Tecplot) in a cuboid (for cut off windows)
int Tecplot_Export::Export_network_threads(const struct cuboid &cub, const vector<vector<Point_3D> > &cnts_points, const int &iwin)const
{
	string file_name = "CNT_Wires_Cutoff_Wins.dat";
	ofstream otec;
	if(iwin==0) 
	{
		otec.open(file_name.c_str());
		otec << "TITLE = CNT_Wires_Cutoff_Wins" << endl;
		otec << "VARIABLES = X, Y, Z" << endl;

		//---------------------------------------------------------------------------
		//Export a 3D cuboid
		if(Export_cuboid(otec, cub)==0) return 0;
	}
	else otec.open(file_name.c_str(), ios_base::app);
	
	//---------------------------------------------------------------------------
	//Export 3D nanotube threads
	if(Export_nano_threads(otec, cub, cnts_points)==0) return 0;

	//---------------------------------------------------------------------------
	otec.close();

	return 1;
}
//---------------------------------------------------------------------------
//Export a 3D cuboid
int Tecplot_Export::Export_cuboid(ofstream &otec, const struct cuboid &cell)const
{
	otec << "ZONE N=" << 8 << ", E=" << 1 << ", F=FEPOINT, ET=BRICK" << endl;
	double bcell_x[2] = {cell.poi_min.x, cell.poi_min.x+cell.len_x};
	double bcell_y[2] = {cell.poi_min.y, cell.poi_min.y+cell.wid_y};
	double bcell_z[2] = {cell.poi_min.z, cell.poi_min.z+0.5*(cell.len_x+cell.wid_y)};
	for(int i=0; i<2; i++)
		for(int j=0; j<2; j++)
			for(int k=0; k<2; k++)
			{
				otec << bcell_x[i] << "  " << bcell_y[j] << "  " << bcell_z[k] << endl;
			}
	
	otec << "1 2 4 3 5 6 8 7" << endl;
	otec << endl << endl;

	return 1;
}
//---------------------------------------------------------------------------
//Export 3D nanotube threads
int Tecplot_Export::Export_nano_threads(ofstream &otec, const struct cuboid &cell, const vector<vector<Point_3D> > &cnts_points)const
{
	//---------------------------------------------------------------------------
	//Export threads
	double cell_x[2] = {cell.poi_min.x, cell.poi_min.x+cell.len_x};
	double cell_y[2] = {cell.poi_min.y, cell.poi_min.y+cell.wid_y};
	double cell_z[2] = {cell.poi_min.z, cell.poi_min.z+cell.hei_z};
	otec << "ZONE N=" << 4 << ", E=" << 1 << ", F=FEPOINT, ET=QUADRILATERAL" << endl;
	for(int i=0; i<2; i++)
		for(int j=0; j<2; j++)
			otec << cell_x[i] << "  " << cell_y[j] << "  " << 0.5*(cell_z[0]+cell_z[1]) << endl;
	otec << "1 2 4 3" << endl;
	otec << endl << endl;

	//---------------------------------------------------------------------------
	//Export threads
	for(int i=0; i<(int)cnts_points.size(); i++)
	{
		otec << "ZONE T=\"Line\"" << endl;
		otec << "i=1," << "j=" << (int)cnts_points[i].size() << ", f=point" << endl;
		for (int j=0; j<(int)cnts_points[i].size(); j++)
		{
			otec << cnts_points[i][j].x << "  " << cnts_points[i][j].y << "  " << cnts_points[i][j].z << endl;
		}
		otec << endl << endl;
	}

	return 1;
}
//---------------------------------------------------------------------------
//The geometric structure of CNT network (by tetrahedron meshes in Tecplot)
int Tecplot_Export::Export_cnt_network_meshes(const struct cuboid &cub, const vector<vector<Point_3D> > &cnts_points, const vector<double> &cnts_radius)const
{
	//For storing the nodes and elements to construct nanotubes with diameters
	vector<vector<Node> > cnts_nodes;
	vector<vector<Element> > cnts_eles;

	//Define a class of GenNetwork for calling for the function below
	GenNetwork *Gentemp = new GenNetwork;
	if(Gentemp->Generate_cnts_nodes_elements(cnts_nodes, cnts_eles, cnts_points, cnts_radius)==0) return 0;
	delete Gentemp;

	//Export nanotube network by tetrahedron elements (Multiple zones in tecplot: each nanotube by one zone)
//	if(Export_cnts_meshes_multizones(cub, cnts_nodes, cnts_eles)==0) return 0;

	//Export nanotube network by tetrahedron elements (Single zones in tecplot: all nanotubes by one zone)
	if(Export_cnts_meshes_singlezone(cub, cnts_nodes, cnts_eles)==0) return 0;

	return 1;
}
//---------------------------------------------------------------------------
//The geometric structure of CNT network for cut off windows (by tetrahedron meshes in Tecplot) 
int Tecplot_Export::Export_cnt_network_meshes(const struct cuboid &cub, const vector<vector<Point_3D> > &cnts_points, const vector<double> &cnts_radius, const int &iwin)const
{
	//For storing the nodes and elements to construct nanotubes with diameters
	vector<vector<Node> > cnts_nodes;
	vector<vector<Element> > cnts_eles;

	//Define a class of GenNetwork for calling for the function below
	GenNetwork *Gentemp = new GenNetwork;
	if(Gentemp->Generate_cnts_nodes_elements(cnts_nodes, cnts_eles, cnts_points, cnts_radius)==0) return 0;
	delete Gentemp;

	//Export nanotube network by tetrahedron elements (Single zones in tecplot: all nanotubes by one zone)
	if(Export_cnts_meshes_singlezone(cub, cnts_nodes, cnts_eles, iwin)==0) return 0;

	return 1;
}
//---------------------------------------------------------------------------
//Export nanotube network by tetrahedron elements (Multiple zones in tecplot: each nanotube by one zone)
int Tecplot_Export::Export_cnts_meshes_multizones(const struct cuboid &cub, const vector<vector<Node> > &nodes, const vector<vector<Element> > &eles)const
{
	//The total number of nanotube threads
	int cnts_account = (int)nodes.size();
	if(cnts_account!=(int)eles.size()) { hout << "Error, the number of node vectors is not the same to the number of element vectors!" << endl; return 0; }
	
	ofstream otec("CNT_Meshes_Multizones.dat");
	otec << "TITLE = CNT_Meshes_Multizones" << endl;
	otec << "VARIABLES = X, Y, Z" << endl;	
	
	//Export a 3D cylinder
	//---------------------------------------------------------------------------
	//Export a 3D cuboid for display in tecplot
	if(Export_cuboid(otec, cub)==0) return 0;

	//---------------------------------------------------------------------------
	//Export a 3D thin cuboid
	otec << "ZONE N=" << 8 << ", E=" << 1 << ", F=FEPOINT, ET=BRICK" << endl;
	double cell_x[2] = {cub.poi_min.x, cub.poi_min.x+cub.len_x};
	double cell_y[2] = {cub.poi_min.y, cub.poi_min.y+cub.wid_y};
	double cell_z[2] = {cub.poi_min.z, cub.poi_min.z+cub.hei_z};
	for(int i=0; i<2; i++)
		for(int j=0; j<2; j++)
			for(int k=0; k<2; k++)
			{
				otec << cell_x[i] << "  " << cell_y[j] << "  " << cell_z[k] << endl;
			}
	
	otec << "1 2 4 3 5 6 8 7" << endl;
	otec << endl << endl;

	//---------------------------------------------------------------------------
	//Export the meshes of nanotubes
	for(int i=0; i<cnts_account; i++)
	{
		otec << "ZONE N=" << (int)nodes[i].size() << ", E=" << (int)eles[i].size() << ", F=FEPOINT, ET=TETRAHEDRON" << endl;
		for(int j=0; j<(int)nodes[i].size(); j++)
		{
			otec << nodes[i][j].x << "  " << nodes[i][j].y << "  " << nodes[i][j].z << endl;
		}
		otec << endl;
		for(int j=0; j<(int)eles[i].size(); j++)
		{
			otec	<< eles[i][j].nodes_id[0]+1 << "  " << eles[i][j].nodes_id[1]+1 << "  " 
					<< eles[i][j].nodes_id[2]+1 << "  " << eles[i][j].nodes_id[3]+1 << endl;
		}
		otec << endl << endl;
	}

	otec.close();

	return 1;
}
//---------------------------------------------------------------------------
//Export nanotube network by tetrahedron elements (Single zones in tecplot: all nanotubes by one zone)
int Tecplot_Export::Export_cnts_meshes_singlezone(const struct cuboid &cub, const vector<vector<Node> > &nodes, const vector<vector<Element> > &eles)const
{
	//The total number of nanotube threads
	int cnts_account = (int)nodes.size();
	if(cnts_account!=(int)eles.size()) { hout << "节点和单元显示的纳米管线数量不一致！ 请检查！" << endl; return 0; }
	
	ofstream otec("CNT_Meshes_Singlezone.dat");
	otec << "TITLE = CNT_Meshes_Singlezone" << endl;
	otec << "VARIABLES = X, Y, Z" << endl;	
	
	//Export a 3D cylinder
	//---------------------------------------------------------------------------
	//Export a 3D cuboid for display in tecplot
	if(Export_cuboid(otec, cub)==0) return 0;

	//---------------------------------------------------------------------------
	//Export a 3D thin cuboid
	otec << "ZONE N=" << 8 << ", E=" << 1 << ", F=FEPOINT, ET=BRICK" << endl;
	double cell_x[2] = {cub.poi_min.x, cub.poi_min.x+cub.len_x};
	double cell_y[2] = {cub.poi_min.y, cub.poi_min.y+cub.wid_y};
	double cell_z[2] = {cub.poi_min.z, cub.poi_min.z+cub.hei_z};
	for(int i=0; i<2; i++)
		for(int j=0; j<2; j++)
			for(int k=0; k<2; k++)
			{
				otec << cell_x[i] << "  " << cell_y[j] << "  " << cell_z[k] << endl;
			}
	
	otec << "1 2 4 3 5 6 8 7" << endl;
	otec << endl << endl;

	//---------------------------------------------------------------------------
	///Export the meshes of nanotubes
	int nodes_num = 0;
	int eles_num = 0;

	for(int i=0; i<cnts_account; i++)
	{
		nodes_num +=  (int)nodes[i].size();
		eles_num += (int)eles[i].size();
	}
		
	otec << "ZONE N=" << nodes_num << ", E=" << eles_num << ", F=FEPOINT, ET=TETRAHEDRON" << endl;
	for(int i=0; i<cnts_account; i++)
	{		
		for(int j=0; j<(int)nodes[i].size(); j++)
		{
			otec << nodes[i][j].x << "  " << nodes[i][j].y << "  " << nodes[i][j].z << endl;
		}
	}
	otec << endl;

	nodes_num = 0;
	for(int i=0; i<cnts_account; i++)
	{
		if(i!=0)  nodes_num +=  (int)nodes[i-1].size();
		for(int j=0; j<(int)eles[i].size(); j++)
		{
			otec	<< eles[i][j].nodes_id[0]+1+nodes_num << "  " << eles[i][j].nodes_id[1]+1+nodes_num << "  " 
					<< eles[i][j].nodes_id[2]+1+nodes_num << "  " << eles[i][j].nodes_id[3]+1+nodes_num << endl;
		}
	}

	otec.close();

	return 1;
}
//---------------------------------------------------------------------------
//Export nanotube network by tetrahedron elements for cut off windows (Single zones in tecplot: all nanotubes by one zone)
int Tecplot_Export::Export_cnts_meshes_singlezone(const struct cuboid &cub, const vector<vector<Node> > &nodes, const vector<vector<Element> > &eles, const int &iwin)const
{
	//The total number of nanotube threads
	int cnts_account = (int)nodes.size();
	if(cnts_account!=(int)eles.size()) { hout << "节点和单元显示的纳米管线数量不一致！ 请检查！" << endl; return 0; }
	
	string file_name = "CNT_Meshes_Singlezone_Cutoff_Wins.dat";
	ofstream otec;
	if(iwin==0) 
	{
		otec.open(file_name.c_str());
		otec << "TITLE = CNT_Meshes_Singlezone_Cutoff_Wins" << endl;
		otec << "VARIABLES = X, Y, Z" << endl;

		//---------------------------------------------------------------------------
		//Export a 3D cuboid
		if(Export_cuboid(otec, cub)==0) return 0;
	}
	else otec.open(file_name.c_str(), ios_base::app);
	
	//---------------------------------------------------------------------------
	//Export a 3D thin cuboid
	otec << "ZONE N=" << 8 << ", E=" << 1 << ", F=FEPOINT, ET=BRICK" << endl;
	double cell_x[2] = {cub.poi_min.x, cub.poi_min.x+cub.len_x};
	double cell_y[2] = {cub.poi_min.y, cub.poi_min.y+cub.wid_y};
	double cell_z[2] = {cub.poi_min.z, cub.poi_min.z+cub.hei_z};
	for(int i=0; i<2; i++)
		for(int j=0; j<2; j++)
			for(int k=0; k<2; k++)
			{
				otec << cell_x[i] << "  " << cell_y[j] << "  " << cell_z[k] << endl;
			}
	
	otec << "1 2 4 3 5 6 8 7" << endl;
	otec << endl << endl;

	//---------------------------------------------------------------------------
	///Export the meshes of nanotubes
	int nodes_num = 0;
	int eles_num = 0;

	for(int i=0; i<cnts_account; i++)
	{
		nodes_num +=  (int)nodes[i].size();
		eles_num += (int)eles[i].size();
	}
		
	otec << "ZONE N=" << nodes_num << ", E=" << eles_num << ", F=FEPOINT, ET=TETRAHEDRON" << endl;
	for(int i=0; i<cnts_account; i++)
	{		
		for(int j=0; j<(int)nodes[i].size(); j++)
		{
			otec << nodes[i][j].x << "  " << nodes[i][j].y << "  " << nodes[i][j].z << endl;
		}
	}
	otec << endl;

	nodes_num = 0;
	for(int i=0; i<cnts_account; i++)
	{
		if(i!=0)  nodes_num +=  (int)nodes[i-1].size();
		for(int j=0; j<(int)eles[i].size(); j++)
		{
			otec	<< eles[i][j].nodes_id[0]+1+nodes_num << "  " << eles[i][j].nodes_id[1]+1+nodes_num << "  " 
					<< eles[i][j].nodes_id[2]+1+nodes_num << "  " << eles[i][j].nodes_id[3]+1+nodes_num << endl;
		}
	}

	otec.close();

	return 1;
}
//---------------------------------------------------------------------------
//The geometric structure of CNT network (by tetrahedron meshes in Tecplot) with a specific filename. This function uses a 1D point vector and a 2D structure vector that references the point vector
int Tecplot_Export::Export_cnt_network_meshes(const struct cuboid &cub, const vector<Point_3D> &cnts_points, const vector<double> &cnts_radius, const vector<vector<long int> > &structure, string filename)const
{
    //For storing the nodes and elements to construct nanotubes with diameters
    vector<vector<Node> > cnts_nodes;
    vector<vector<Element> > cnts_eles;
    
    //Define a class of GenNetwork for calling for the function below
    GenNetwork *Gentemp = new GenNetwork;
    if(Gentemp->Generate_cnts_nodes_elements(cnts_nodes, cnts_eles, cnts_points, cnts_radius, structure)==0) return 0;
    delete Gentemp;
    
    //Export nanotube network by tetrahedron elements (Multiple zones in tecplot: each nanotube by one zone)
    if(Export_cnts_meshes_singlezone(cub, cnts_nodes, cnts_eles, filename)==0) return 0;
    
    return 1;
}
//---------------------------------------------------------------------------
//Export nanotube network by tetrahedron elements (Multiple zones in tecplot: each nanotube by one zone) with a specific filename
int Tecplot_Export::Export_cnts_meshes_multizones(const struct cuboid &cub, const vector<vector<Node> > &nodes, const vector<vector<Element> > &eles, string filename)const
{
    //The total number of nanotube threads
    int cnts_account = (int)nodes.size();
    if(cnts_account!=(int)eles.size()) { hout << "Error, the number of node vectors is not the same to the number of element vectors!" << endl; return 0; }
    
    ofstream otec(filename.c_str());
    otec << "TITLE = CNT_Meshes_Multizones" << endl;
    otec << "VARIABLES = X, Y, Z" << endl;
    
    //---------------------------------------------------------------------------
    //Export a 3D cylinder
    if(Export_cuboid(otec, cub)==0) return 0;
    
    //---------------------------------------------------------------------------
    //Export the meshes of nanotubes
    for(int i=0; i<cnts_account; i++)
    {
        otec << "ZONE N=" << (int)nodes[i].size() << ", E=" << (int)eles[i].size() << ", F=FEPOINT, ET=TETRAHEDRON" << endl;
        for(int j=0; j<(int)nodes[i].size(); j++)
        {
            otec << nodes[i][j].x << "  " << nodes[i][j].y << "  " << nodes[i][j].z << endl;
        }
        otec << endl;
        for(int j=0; j<(int)eles[i].size(); j++)
        {
            otec	<< eles[i][j].nodes_id[0]+1 << "  " << eles[i][j].nodes_id[1]+1 << "  "
            << eles[i][j].nodes_id[2]+1 << "  " << eles[i][j].nodes_id[3]+1 << endl;
        }
        otec << endl << endl;
    }
    
    otec.close();
    
    return 1;
}
//---------------------------------------------------------------------------
//Export nanotube network by tetrahedron elements (Single zones in tecplot: all nanotubes by one zone) with a specific filename
int Tecplot_Export::Export_cnts_meshes_singlezone(const struct cuboid &cub, const vector<vector<Node> > &nodes, const vector<vector<Element> > &eles, string filename)const
{
    //The total number of nanotube threads
    int cnts_account = (int)nodes.size();
    if(cnts_account!=(int)eles.size()) { hout << "节点和单元显示的纳米管线数量不一致！ 请检查！" << endl; return 0; }
    
    ofstream otec(filename.c_str());
    otec << "TITLE = CNT_Meshes_Singlezone" << endl;
    otec << "VARIABLES = X, Y, Z" << endl;
    
    //---------------------------------------------------------------------------
    //Export a 3D cylinder
    if(Export_cuboid(otec, cub)==0) return 0;
    
    //---------------------------------------------------------------------------
    ///Export the meshes of nanotubes
    int nodes_num = 0;
    int eles_num = 0;
    
    for(int i=0; i<cnts_account; i++)
    {
        nodes_num +=  (int)nodes[i].size();
        eles_num += (int)eles[i].size();
    }
    
    otec << "ZONE N=" << nodes_num << ", E=" << eles_num << ", F=FEPOINT, ET=TETRAHEDRON" << endl;
    for(int i=0; i<cnts_account; i++)
    {
        for(int j=0; j<(int)nodes[i].size(); j++)
        {
            otec << nodes[i][j].x << "  " << nodes[i][j].y << "  " << nodes[i][j].z << endl;
        }
    }
    otec << endl;
    
    nodes_num = 0;
    for(int i=0; i<cnts_account; i++)
    {
        if(i!=0)  nodes_num +=  (int)nodes[i-1].size();
        for(int j=0; j<(int)eles[i].size(); j++)
        {
            otec	<< eles[i][j].nodes_id[0]+1+nodes_num << "  " << eles[i][j].nodes_id[1]+1+nodes_num << "  "
            << eles[i][j].nodes_id[2]+1+nodes_num << "  " << eles[i][j].nodes_id[3]+1+nodes_num << endl;
        }
    }
    
    otec.close();
    
    return 1;
}
//===========================================================================
