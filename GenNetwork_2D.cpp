//====================================================================================
//SOFTWARE:	Network of Eelectrically Conductive Nanocomposites (NECN)
//CODE FILE:	GenNetwork.cpp
//OBJECTIVE:	To generate networks with overlapping
//AUTHOR:		Fei Han; Angel Mora
//E-MAIL:			fei.han@kaust.edu.sa	;	angel.mora@kaust.edu.sa
//====================================================================================

#include "GenNetwork_2D.h"
#include "Geometry_2D.h"
#include "Geometry_3D.h"

//Generate 2D nanowire networks with ovelapping
int GenNetwork::Generate_nanowire_networks(const struct Geom_RVE &geom_rve, const struct Nanowire_Geo &nanowire_geo, vector<Point_2D> &cpoints, vector<double> &cnts_radius, vector<vector<int> > &cstructures)const
{
	// Define a vector of vectors for storing the cartesian coordinates of nanowires
    vector<vector<Point_2D> > cnts_points;
	//Use the Mersenne Twister for the random number generation
	if (Generate_network_threads_mt(geom_rve, nanowire_geo, cnts_points, cnts_radius)==0) return 0;
    //-----------------------------------------------------------------------------------------------------------------------------------------
    //Transform the 2D cnts_points into 1D cpoints and 2D cstructuers
    if (Transform_cnts_points(cnts_points, cpoints, cstructures)==0) return 0;	
    //-----------------------------------------------------------------------------------------------------------------------------------------
    //A new class of Tecplot_Export
    Tecplot_Export *Tecexpt = new Tecplot_Export;
     
	//Generate a rectangle for RVE
	Rectangle exrect(geom_rve.ex_origin, geom_rve.ex_len, geom_rve.ey_wid);
     
     //The geometric structure of CNT network (by threads in Tecplot)
     if(Tecexpt->Export_network_threads_2D(exrect, cnts_points)==0) return 0;
     
	 //Find the max cnt radius
	 double max_rad = 0;
	 for(int i=0; i<(int)cnts_radius.size(); i++)
	 {
		 if(max_rad<cnts_radius[i]) max_rad = cnts_radius[i];
	 }

	 //Generate a cuboid for RVE
	 struct cuboid cub;
     cub.poi_min.x = geom_rve.ex_origin.x;
     cub.poi_min.y = geom_rve.ex_origin.y;
	 cub.poi_min.z = 0.0;
     cub.len_x = geom_rve.ex_len;
     cub.wid_y = geom_rve.ey_wid;
	 cub.hei_z = 2*max_rad;

	//Define a vector of vectors for qusi2D output
    vector<vector<Point_3D> > cnts_3Dpois;
	 for(int i=0; i<(int)cnts_points.size(); i++)
	 {
		vector<Point_3D> vec_3dps;
		//first point
		Point_3D temp_poi(cnts_points[i][0].x, cnts_points[i][0].y, max_rad);
		vec_3dps.push_back(temp_poi);

		//middle points
		//double x0 = cnts_points[i][0].x;
		//double deltx = cnts_points[i][1].x - cnts_points[i][0].x;
		//double y0 = cnts_points[i][0].y;
		//double delty = cnts_points[i][1].y - cnts_points[i][0].y;
		//for(int j=1; j<10; j++)
		//{
		//	temp_poi.x = x0+deltx*j*0.1;
		//	temp_poi.y = y0+delty*j*0.1;
		//	temp_poi.z = max_rad;
		//	vec_3dps.push_back(temp_poi);
		//}

		//end point
		temp_poi.x = cnts_points[i][1].x;
		temp_poi.y = cnts_points[i][1].y;
		temp_poi.z = max_rad;
		vec_3dps.push_back(temp_poi);
		 
		//push_back
		cnts_3Dpois.push_back(vec_3dps);
	 }

     //The geometric structure of CNT network (by tetrahedron meshes in Tecplot)
	 //Attention: little parts of nanotube volumes out of the rectangle
     if(Tecexpt->Export_cnt_network_meshes(cub, cnts_3Dpois, cnts_radius)==0) return 0;

	return 1;
}	
//Generate a network defined by points and connections
//Use the Mersenne Twister for the random number generation
int GenNetwork::Generate_network_threads_mt(const struct Geom_RVE &geom_rve, const struct Nanowire_Geo &nanowire_geo, vector<vector<Point_2D> > &Cnts_points,  vector<double> &Cnts_radius)const
{
    //Generate random seed in terms of local time
    //unsigned int time_seed = 1453384844;
    //unsigned int time_seed = (unsigned int)time(NULL);
    //hout << "Time seed "<<time_seed<<endl;
    //srand(time_seed);
    //srand((unsigned int)time(NULL));
    
    //---------------------------------------------------------------------------
    //Set up the Mersenne Twisters used for the different variables
    //Use random_device to generate a seed for Mersenne twister engine.
    std::random_device rd;
    // Use Mersenne twister engine to generate pseudo-random numbers.
    //Generate differnet engines for different variables
    std::mt19937 engine_x(rd());
    std::mt19937 engine_y(rd());
    std::mt19937 engine_pha(rd());
	
    //"Filter" MT's output to generate double values, uniformly distributed on the closed interval [0, 1].
    std::uniform_real_distribution<double> dist(0.0, 1.0);
	
    double area_sum = 0;  //the sum of area of generated CNTs
    double wt_sum = 0;   //the sum of weight of generated CNTs
    int cnt_seed_count = 0; //to record the number of generated seeds of a CNT 
    int point_overlap_count = 0; //to record the number of overlappings to calculate the projected area properly. This gives the transparency %T=100-a1*A_nw
    int point_overlap_count_unique = 0; //to record the number of points that were overlapping other points ******May not be necessary, not sure : we will see later
    const int MAX_ATTEMPTS = 5; // Overlapping case fixing
	
	// ************** I suppose the following vectors are used to handle the overlappings *******************
	// Global coordinates global_coordinates[0].at(i) -> global_coordinates[1].at(i) i'th nanowire
    vector<vector<int> > global_coordinates;
    //sectioned_domain[i] contains all the points in sub-region i.
    //Sub-region i is an overlapping subregion to check for penetrations
    vector<vector<long int> > sectioned_domain;
    //n_subregions[0] is the number of subregions along x
    //n_subregions[1] is the number of subregions along y
    vector<int> n_subregions;
    //Initialize the vector sub-regions
    Initialize_subregions(geom_rve, n_subregions, sectioned_domain);
	
	
	
    //Generate Rectangles that represent the extended domain and the composite domain
    //To calculate the effective portion (length) which falls into the given region (RVE)
	//generate a Rectangle to represent the composite domain
	Rectangle gvcub(geom_rve.origin, geom_rve.len_x, geom_rve.wid_y);
	//generate a Rectangle to represent the extended domain
	Rectangle excub(geom_rve.ex_origin, geom_rve.ex_len, geom_rve.ey_wid);

	//--------------- GENERATING SEEDS RANDOMLY AND THEN THE END POINT UNTIL AREA_SUM IS REACHED -----------
    while((nanowire_geo.criterion == "area"&&area_sum < nanowire_geo.real_area)||
			(nanowire_geo.criterion == "wt"&&wt_sum < nanowire_geo.real_weight))
    {
        //---------------------------------------------------------------------------
        //Define two points for a new nanowire
        Point_2D new_cnt[2];
        
        //---------------------------------------------------------------------------
        //Randomly generate a length of a CNT
        double cnt_length;
        if(Get_random_value_mt(nanowire_geo.len_distrib_type, engine_pha, dist, nanowire_geo.len_min, nanowire_geo.len_max, cnt_length)==0) return 0;
        //---------------------------------------------------------------------------

        //Randomly generate a radius of a CNT
        double cnt_rad;
        if(Get_random_value_mt(nanowire_geo.rad_distrib_type, engine_pha, dist, nanowire_geo.rad_min, nanowire_geo.rad_max, cnt_rad)==0) return 0;
        
        //---------------------------------------------------------------------------
        //Randomly generate the orientation of the nanowire 'Pha' is the orintation of the nanowire (-PI/2, PI/2)
        double cnt_pha;
        if(Get_uniform_direction_mt(nanowire_geo, cnt_pha, engine_pha, dist)==0) return 0;
		
        //---------------------------------------------------------------------------
        //The increased weight of each nanowire (If the different radii of nanowire are considered, the linear_density may be different in every nanowire)
        const double wei_para = nanowire_geo.linear_density;
		
        //---------------------------------------------------------------------------
        //Randomly generate a seed (initial point) of a CNT in the extended RVE	
        
        Point_2D cnt_poi;
        if(Get_seed_point_mt(excub, cnt_poi, engine_x, engine_y, dist)==0) return 0;
		
		//int counter=1;
  //      //Check overlapping of the initial point 	
		//while(!Check_overlapping(Cnts_points, cnt_poi)) {
		//	if(Get_seed_point_mt(excub, cnt_poi, engine_x, engine_y, dist)==0) return 0;
  //          cnt_seed_count++;					//record the number of seed generations
  //          //hout << "Seed deleted" << endl;
  //          if (counter == MAX_ATTEMPTS) {
  //              hout << "Too many attempts to resolve overlapping of an intial CNT point (" << counter << " attempts). ";
  //              hout << cnt_poi.x << ' ' << cnt_poi.y << endl;
  //              return 0;
  //          }
  //          counter ++;
		//} 
		
        new_cnt[0] = cnt_poi;	//store this seed point in the array for a new nanowire
        //---------------------------------------------------------------------------
        cnt_seed_count++;					//record the number of seed generations
        double max_seed = 1E13;
        if(cnt_seed_count>max_seed)
        {
            hout << "The number of seed genrations is lager than "<<max_seed<<", but the nanowire generation still fails to acheive the demanded area fraction." << endl;
            return 0;
        }
		
		// We have generated a new seed (Initial point of the nanowire) and the orientation of the nanowire so, the end point of nanowire is
		Point_2D point_add(cnt_length*cos(cnt_pha),cnt_length*sin(cnt_pha));
		cnt_poi = cnt_poi + point_add;
		// We will save this as the end point only if it doesnt go out of the RVE 
		
        //---------------------------------------------------------------------------
        //If the new CNT point grows out of the RVE, the intersecting point at the edge of RVE will be calculated.
        //The length going out of RVE will be cut and (the outside part will be translated into the RVE for the periodical boundary condition) @Fei How is it translated? I don't see it in the code           
		if(excub.contain_in(cnt_poi)==0)
        {
			Point_2D touch_edge;  //intersection point
			if(Get_intersecting_point_RVE_edge(excub, new_cnt[0], cnt_poi, touch_edge)==0) 
			{
				hout << "Error, it fails to find the intersecting point of RVE edge."<<endl;
				return 0;
			}
			cnt_poi = touch_edge;
		}

		//---------------------------------------------------------------------------
        // At this point a nanowire is created with length 'cnt_length' and diameter '2*cnt_rad' and orientation 'cnt_pha'
	    // Calculate the accumulated area and the weight 
        double temp_length = Effective_length_given_region(gvcub, new_cnt[0],cnt_poi);
        if (temp_length > 0.0)
        {
            area_sum += temp_length*2*cnt_rad;		//an accumulation on the area (Not the projection area)
            wt_sum += temp_length*wei_para;	     	//an accumulation on the weight
        }

        //---------------------------------------------------------------------------
	    // Add the point as the end point of the nanowire
	    new_cnt[1] = cnt_poi;
		        
        //---------------------------------------------------------------------------
        //Store the CNT points
		//If the new_cnt array has exactly two points - error if something else happens
		vector<Point_2D> cnt_temp;
		cnt_temp.push_back(new_cnt[0]);  // First point of the nanowire
		cnt_temp.push_back(new_cnt[1]);  // End point of the same nanowire
		Cnts_points.push_back(cnt_temp);
        Cnts_radius.push_back(cnt_rad);
		cnt_temp.clear();
	}	
    if(nanowire_geo.criterion == "area") hout << "    The area fraction of generated CNTs is about : " << area_sum/geom_rve.area << endl;
    
//  hout << "There were " << point_overlap_count_unique << " overlapping points and ";
//  hout << point_overlap_count << " overlaps, " << endl;
    
    return 1;
}

//---------------------------------------------------------------------------
//Generate a random value through a probability distribution function
int GenNetwork::Get_random_value_mt(const string &dist_type, mt19937 &engine, uniform_real_distribution<double> &dist, const double &min, const double &max, double &value)const
{
    //Check if limits are correctly defined
    if(min>max) { hout << "Error, the minimum value is larger than the maximum value (Get_random_value)!" << endl; return 0; }
    
    //Check if the interval has 0 length
    if (max == min) {
        //In this case, value is either of the limits
        //To be consistent with the formulation below, value is set equal to min
        value = min;
        return 1;
    }
    
    if(dist_type=="uniform")	//uniform distribution
    {
        value = (max-min)*dist(engine) + min;
    }
    else if(dist_type=="normal")	//normal distribution
    {
        double sum=0;
        for(int i=0; i<12; i++)
        {
            sum = sum + dist(engine);
        }
        value = (max-min)*sum/12.0 + min;
    }
    
    return 1;
}

//---------------------------------------------------------------------------
//Randomly generate a seed (intial point) of a CNT in the RVE
int GenNetwork::Get_seed_point_mt(const Rectangle &cub, Point_2D &poin, mt19937 &engine_x, mt19937 &engine_y, uniform_real_distribution<double> &dist)const
{
    
    poin.x = cub.point[0].x + cub.length*dist(engine_x);
    
    poin.y = cub.point[0].y+ cub.width*dist(engine_y);
    
    return 1;
}

//---------------------------------------------------------------------------
//Randomly generate a direction for a given nanowire 
int GenNetwork::Get_uniform_direction_mt(const struct Nanowire_Geo &nanowire_geo, double &cnt_pha, mt19937 &engine_pha, uniform_real_distribution<double> &dist)const
{
    if(nanowire_geo.dir_distrib_type=="random") //random distribution
    {   
        cnt_pha = nanowire_geo.angle_min + (nanowire_geo.angle_max-nanowire_geo.angle_min) * dist(engine_pha);
		if(dist(engine_pha)<0.5) cnt_pha = -cnt_pha;
    }
    else if(nanowire_geo.dir_distrib_type=="normal") //normal distribution
    {
        double sum=0;
        for(int i=0; i<12; i++)
        {
            sum = sum + dist(engine_pha);
        }
		cnt_pha = (nanowire_geo.angle_max-nanowire_geo.angle_min)*sum/12.0 + nanowire_geo.angle_min;
		if(dist(engine_pha)<0.5) cnt_pha = -cnt_pha;
	}

	return 1;
}

//---------------------------------------------------------------------------
double GenNetwork::Effective_length_given_region(const Rectangle &cub, const Point_2D &last_point, const Point_2D &new_point)const
{
    //Check if the last point is inside the given region
    int last_bool = cub.contain_in(last_point);
    //Check if the new point is inside the given region
    int new_bool = cub.contain_in(new_point);
    
    //To store the intersecting point
    Point_2D touchpoint;
    
    //Decide the corresponding case and calculate area fraction
    if (last_bool&&new_bool)
        return last_point.distance_to(new_point); //both points are inside so add the total length
    else if (last_bool&&(!new_bool))  //if the last point is inside and the new point is outside
    {
        if(Get_intersecting_point_RVE_edge(cub, last_point, new_point, touchpoint)==0){
            hout << "Error in Effective_length_given_region, case last_bool&&(!new_bool) "<<endl;
            return 0;
        }
        return last_point.distance_to(touchpoint);
    }
    else if ((!last_bool)&&new_bool)  //if the last point is outside and the new point is inside
    {
        if(Get_intersecting_point_RVE_edge(cub, new_point, last_point, touchpoint)==0) {
            hout << "Error in Effective_length_given_region, case (!last_bool)&&new_bool"<<endl;
            return 0;
        }
        return new_point.distance_to(touchpoint);
    }
    else
        return 0.0; //if both points are outside
}

//---------------------------------------------------------------------------
//This functions initializes the vectors n_subregions and sectioned_domain
//
//The n_subregions vector is defined to avoid calculating the number of sub-regions for every point when the functions
//Default_region and Add_to_subregion are called. Thus, saving computational time.
//n_subregions[0] is the number of subregions along x
//n_subregions[1] is the number of subregions along y

//The vector sectioned_domain contains the sub-regions to look for overlapping
//It is initialized with the number of sub-regions in the sample
void GenNetwork::Initialize_subregions(const struct Geom_RVE &geom_rve, vector<int> &nsubregions, vector<vector<long int> > &sectioned_domain)const
{
    //Initialize nsubregions
    
    //variable to store the number of subregions
    int s;
    //Number of subregions along x
    s = (int)(geom_rve.len_x/geom_rve.gs_minx);
    //Add the number of sub-regions to the vector
    nsubregions.push_back(s);
    //Number of subregions along y
    s = (int)(geom_rve.wid_y/geom_rve.gs_miny);
    //Add the number of sub-regions to the vector
    nsubregions.push_back(s);
    
    //Initialize sectioned_domain
    vector<long int> empty;
    sectioned_domain.assign(nsubregions[0]*nsubregions[1], empty);
}

//---------------------------------------------------------------------------
//This function returns the subregion a point belongs to
int GenNetwork::Get_subregion(const struct Geom_RVE &geom_rve, const vector<int> &n_subregions, const Point_2D &point)const
{
    if (Judge_RVE_including_point(geom_rve, point)) {
        //These variables will give me the region cordinates of the region that a point belongs to
        int a, b;
        //Calculate the region-coordinates
        a = (int)((point.x-geom_rve.origin.x)/geom_rve.gs_minx);
        //Limit the value of a as it has to go from 0 to n_subregions[0]-1
        if (a == n_subregions[0]) a--;
        b = (int)((point.y-geom_rve.origin.y)/geom_rve.gs_miny);
        //Limit the value of b as it has to go from 0 to n_subregions[1]-1
        if (b == n_subregions[1]) b--;
        return (a + (b*n_subregions[0]) );
    } else {
        //If the point is in the boundary layer, then there is no need to calculate sub-region
        return -1;
    }
}

//---------------------------------------------------------------------------
//To judge if a point is included in a RVE
int GenNetwork::Judge_RVE_including_point(const struct Geom_RVE &geom_rve, const Point_2D &point)const
{
    if(point.x<geom_rve.origin.x||point.x>geom_rve.origin.x+geom_rve.len_x||
       point.y<geom_rve.origin.y||point.y>geom_rve.origin.y+geom_rve.wid_y) 
	   return 0;
    
    return 1;
}

//---------------------------------------------------------------------------
 int GenNetwork::Check_overlapping(const vector<vector<Point_2D> > &cnts_points, Point_2D cnt_poi)const		
{
	 for (unsigned int i = 0; i < cnts_points.size(); i++) {
		 for (int j = 0; j < (int)cnts_points[i].size(); j++) {
			 if (cnt_poi == cnts_points[i][j] ) {
				 // New generated seed overlaps with one of the nanowires end points already generated
				 return -1;
			 }
		 }
	 }
	return 1;
}

 //---------------------------------------------------------------------------
int GenNetwork::Get_intersecting_point_RVE_edge(const Rectangle &cub, const Point_2D &point1, const Point_2D &point2, Point_2D &touch_point)const
{
	// You have 4 edges edge1: left vertical line other edges-- Anti clockwise
	// Starting from edge 1, we will check for intersection one by one
	// for example, Checking for intersection with edge 2 defined by (Rectangle.point[0],Rectangle.point[3])
	// AB is a line formed by p0p1; CD is a line formed by p2p3; we are finding the intersection of segments AB and CD
	double p0_x, p0_y, p1_x, p1_y, p2_x, p2_y, p3_x, p3_y;
	
	p0_x = point1.x;   p0_y = point1.y;  // A
	p1_x = point2.x;   p1_y = point2.y;  // B
	
		int i1, i2; // Two points of an edge
	for (int j=0; j<4; j++)	{
		if (j==0) {i1=1; i2=0;}  // Left edge
		else if (j==1) {i1=0; i2=3;}  // Bottom edge
		else if (j==2) {i1=3; i2=2;}  // right edge
		else if (j==3) {i1=2; i2=1;}  // Top edge 
		p2_x = cub.point[i1].x; 	p2_y = cub.point[i1].y;  // C
		p3_x = cub.point[i2].x; 	p3_y = cub.point[i2].y;  // D
		if (get_line_intersection(p0_x, p0_y, p1_x, p1_y, p2_x, p2_y, p3_x, p3_y, touch_point) ==1) {
		// Intersection found, now stop the loop 
		return 1;
		}
		}
	return 0; 
}	

//---------------------------------------------------------------------------
int GenNetwork::get_line_intersection(const double &p0_x, const double &p0_y, const double &p1_x, const double &p1_y, const double &p2_x, const double &p2_y, const double &p3_x, const double &p3_y, Point_2D &touchpoint)const
{
    double s1_x, s1_y, s2_x, s2_y;
    s1_x = p1_x - p0_x;     s1_y = p1_y - p0_y;
    s2_x = p3_x - p2_x;     s2_y = p3_y - p2_y;

    double s, t;
    s = (-s1_y * (p0_x - p2_x) + s1_x * (p0_y - p2_y)) / (-s2_x * s1_y + s1_x * s2_y);
    t = ( s2_x * (p0_y - p2_y) - s2_y * (p0_x - p2_x)) / (-s2_x * s1_y + s1_x * s2_y);

    if (s >= 0 && s <= 1 && t >= 0 && t <= 1)
    {
        // Collision detected
        if (touchpoint.x != NULL)
            touchpoint.x = p0_x + (t * s1_x);
        if (touchpoint.y != NULL)
            touchpoint.y = p0_y + (t * s1_y);
        return 1;
    }
    return 0; // No collision
}

//---------------------------------------------------------------------------
//Transform the 2D cnts_points into 1D cpoints and 2D cstructuers
int GenNetwork::Transform_cnts_points(const vector<vector<Point_2D> > &cnts_points, vector<Point_2D> &cpoints, vector<vector<int> > &cstructures)const
{
    hout << "There are "<<cnts_points.size()<<" CNTs."<<endl;
    long int count = 0;
    for(int i=0; i<(int)cnts_points.size(); i++)
    {
        vector<int> struct_temp;
        for(int j=0; j<(int)cnts_points[i].size(); j++) // right now we are only considering the end points of the nanowires, later we will split the NWs
        {
            cpoints.push_back(cnts_points[i][j]);
            cpoints.back().flag = i;
            struct_temp.push_back(count);
            count++;
        }
        cstructures.push_back(struct_temp);
    }
    
    hout << "There are "<<cpoints.size() << " points."<<endl;
    
    return 1;
}

//---------------------------------------------------------------------------
//Generate the nodes and tetrahedron elements of nanotubes (No const following this function because a sum operation on two Point_3D points inside)
int GenNetwork::Generate_cnts_nodes_elements(vector<vector<Node> > &nodes, vector<vector<Element> > &eles, const vector<vector<Point_3D> > &cnts_points, const vector<double> &cnts_radius)
{
    //Looping the generated nanotubes
    for(int i=0; i<(int)cnts_points.size(); i++)
    {
        vector<Node> nod_temp;
        vector<Element> ele_temp;
        
        const int cps = (int)cnts_points[i].size();
        for(int j=0; j<cps; j++)
        {
            //Calculate the normal vector of the plane
            Point_3D plane_normal;
            if(j==0)
            {
                plane_normal.x = cnts_points[i][j].x - cnts_points[i][j+1].x;
                plane_normal.y = cnts_points[i][j].y - cnts_points[i][j+1].y;
                plane_normal.z = cnts_points[i][j].z - cnts_points[i][j+1].z;
            }
            else if(j==cps-1)
            {
                plane_normal.x = cnts_points[i][j-1].x - cnts_points[i][j].x;
                plane_normal.y = cnts_points[i][j-1].y - cnts_points[i][j].y;
                plane_normal.z = cnts_points[i][j-1].z - cnts_points[i][j].z;
            }
            else
            {
                const Point_3D vect[3] = { cnts_points[i][j-1], cnts_points[i][j], cnts_points[i][j+1] };
                const double A = pow(vect[0].x-vect[1].x,2)+pow(vect[0].y-vect[1].y,2)+pow(vect[0].z-vect[1].z,2);
                const double B = pow(vect[2].x-vect[1].x,2)+pow(vect[2].y-vect[1].y,2)+pow(vect[2].z-vect[1].z,2);
                const double tt = sqrt(A/B);
                
                //Coordinate transformation to find the intersection points
                double x, y, z;
                x=vect[1].x+tt*(vect[2].x-vect[1].x);
                y=vect[1].y+tt*(vect[2].y-vect[1].y);
                z=vect[1].z+tt*(vect[2].z-vect[1].z);
                
                plane_normal.x = vect[0].x - x;
                plane_normal.y = vect[0].y - y;
                plane_normal.z = vect[0].z - z;
            }
            //The center point of the circle on the plane
            Point_3D plane_center = cnts_points[i][j];
            
            //Define the number of sections along the circumference
            const int num_sec = 36;
            if(j==0)
            {
                double normal_sita, normal_pha;  //Direction angles
                //Calculate the angles of the normal verctor of the plane in the spherical coordinate
                if(Get_angles_vector_in_spherial_coordinates(plane_normal, normal_sita, normal_pha)==0) return 0;
                
                //Calculate a group of equidistant points along the circumference which is on the plane defined by the center point of the circle and the normal vector
                if(Get_points_circle_in_plane(plane_center, normal_sita, normal_pha, cnts_radius[i], num_sec, nod_temp)==0) return 0;
            }
            else
            {
                //Calculate a group of projected points (which are on the plane with the center point of the circle and the normal vector)
                //which are projected from a group of points on the previous circumference and projected along the direction of line_vec
                Point_3D line_vec;
                line_vec.x = cnts_points[i][j-1].x - cnts_points[i][j].x;
                line_vec.y = cnts_points[i][j-1].y - cnts_points[i][j].y;
                line_vec.z = cnts_points[i][j-1].z - cnts_points[i][j].z;
                if(Get_projected_points_in_plane(plane_center, plane_normal, line_vec, num_sec, nod_temp)==0) return 0;
            }
            
            //Generate a vector of elements
            if(j!=0)
            {
                int nodes_num[6];
                nodes_num[0] = (j-1)*(num_sec+1);   //The number of the center
                nodes_num[3] = j*(num_sec+1);
                for(int k=1; k<=num_sec; k++)
                {
                    nodes_num[1] = (j-1)*(num_sec+1) + k;
                    nodes_num[2] = (j-1)*(num_sec+1) + 1 + k%num_sec;
                    nodes_num[4] = j*(num_sec+1) + k;
                    nodes_num[5] = j*(num_sec+1) + 1 + k%num_sec;
                    
                    Element eles_num[3];
                    //----------------------------------------------------------------
                    //Insert the numbers of nodes to the elements
                    eles_num[0].nodes_id.push_back(nodes_num[0]);
                    eles_num[0].nodes_id.push_back(nodes_num[1]);
                    eles_num[0].nodes_id.push_back(nodes_num[2]);
                    eles_num[0].nodes_id.push_back(nodes_num[3]);
                    
                    eles_num[1].nodes_id.push_back(nodes_num[1]);
                    eles_num[1].nodes_id.push_back(nodes_num[2]);
                    eles_num[1].nodes_id.push_back(nodes_num[3]);
                    eles_num[1].nodes_id.push_back(nodes_num[5]);
                    
                    eles_num[2].nodes_id.push_back(nodes_num[1]);
                    eles_num[2].nodes_id.push_back(nodes_num[3]);
                    eles_num[2].nodes_id.push_back(nodes_num[4]);
                    eles_num[2].nodes_id.push_back(nodes_num[5]);
                    //----------------------------------------------------------------
                    //Insert the number of elements to element vector
                    ele_temp.push_back(eles_num[0]);
                    ele_temp.push_back(eles_num[1]);
                    ele_temp.push_back(eles_num[2]);
                }
            }
        }
        
        nodes.push_back(nod_temp);
        eles.push_back(ele_temp);
    }
    
    return 1;
}

//---------------------------------------------------------------------------
//Generate the nodes and tetrahedron elements of nanotubes (No const following this function because a sum operation on two Point_3D points inside). This function uses a 1D point vector and a 2D structure vector that references the point vector
int GenNetwork::Generate_cnts_nodes_elements(vector<vector<Node> > &nodes, vector<vector<Element> > &eles, const vector<Point_3D> &cnts_points, const vector<double> &cnts_radius, const vector<vector<long int> > &structure)
{
    //Looping the generated nanotubes
    for(int i=0; i<(int)structure.size(); i++)
    {
        vector<Node> nod_temp;
        vector<Element> ele_temp;
        
        const int cps = (int)structure[i].size();
        for(int j=0; j<cps; j++)
        {
            //Calculate the normal vector of the plane
            Point_3D plane_normal;
            if(j==0)
            {
                long int P1 = structure[i][j];
                long int P2 = structure[i][j+1];
                plane_normal.x = cnts_points[P1].x - cnts_points[P2].x;
                plane_normal.y = cnts_points[P1].y - cnts_points[P2].y;
                plane_normal.z = cnts_points[P1].z - cnts_points[P2].z;
            }
            else if(j==cps-1)
            {
                long int P1 = structure[i][j-1];
                long int P2 = structure[i][j];
                plane_normal.x = cnts_points[P1].x - cnts_points[P2].x;
                plane_normal.y = cnts_points[P1].y - cnts_points[P2].y;
                plane_normal.z = cnts_points[P1].z - cnts_points[P2].z;
            }
            else
            {
                long int P1 = structure[i][j-1];
                long int P2 = structure[i][j];
                long int P3 = structure[i][j+1];
                const Point_3D vect[3] = { cnts_points[P1], cnts_points[P2], cnts_points[P3] };
                const double A = pow(vect[0].x-vect[1].x,2)+pow(vect[0].y-vect[1].y,2)+pow(vect[0].z-vect[1].z,2);
                const double B = pow(vect[2].x-vect[1].x,2)+pow(vect[2].y-vect[1].y,2)+pow(vect[2].z-vect[1].z,2);
                const double tt = sqrt(A/B);
                
                //Coordinate transformation to find the intersection points
                double x, y, z;
                x=vect[1].x+tt*(vect[2].x-vect[1].x);
                y=vect[1].y+tt*(vect[2].y-vect[1].y);
                z=vect[1].z+tt*(vect[2].z-vect[1].z);
                
                plane_normal.x = vect[0].x - x;
                plane_normal.y = vect[0].y - y;
                plane_normal.z = vect[0].z - z;
            }
            //The center point of the circle on the plane
            Point_3D plane_center = cnts_points[structure[i][j]];
            
            //Define the number of sections along the circumference
            const int num_sec = 36;
            if(j==0)
            {
                double normal_sita, normal_pha;  //Direction angles
                //Calculate the angles of the normal verctor of the plane in the spherical coordinate
                if(Get_angles_vector_in_spherial_coordinates(plane_normal, normal_sita, normal_pha)==0) return 0;
                
                //Calculate a group of equidistant points along the circumference which is on the plane defined by the center point of the circle and the normal vector
                if(Get_points_circle_in_plane(plane_center, normal_sita, normal_pha, cnts_radius[i], num_sec, nod_temp)==0) return 0;
            }
            else
            {
                //Calculate a group of projected points (which are on the plane with the center point of the circle and the normal vector)
                //which are projected from a group of points on the previous circumference and projected along the direction of line_vec
                Point_3D line_vec;
                long int P1 = structure[i][j-1];
                long int P2 = structure[i][j];
                line_vec.x = cnts_points[P1].x - cnts_points[P2].x;
                line_vec.y = cnts_points[P1].y - cnts_points[P2].y;
                line_vec.z = cnts_points[P1].z - cnts_points[P2].z;
                if(Get_projected_points_in_plane(plane_center, plane_normal, line_vec, num_sec, nod_temp)==0) return 0;
            }
            
            //Generate a vector of elements
            if(j!=0)
            {
                int nodes_num[6];
                nodes_num[0] = (j-1)*(num_sec+1);   //The number of the center
                nodes_num[3] = j*(num_sec+1);
                for(int k=1; k<=num_sec; k++)
                {
                    nodes_num[1] = (j-1)*(num_sec+1) + k;
                    nodes_num[2] = (j-1)*(num_sec+1) + 1 + k%num_sec;
                    nodes_num[4] = j*(num_sec+1) + k;
                    nodes_num[5] = j*(num_sec+1) + 1 + k%num_sec;
                    
                    Element eles_num[3];
                    //----------------------------------------------------------------
                    //Insert the numbers of nodes to the elements
                    eles_num[0].nodes_id.push_back(nodes_num[0]);
                    eles_num[0].nodes_id.push_back(nodes_num[1]);
                    eles_num[0].nodes_id.push_back(nodes_num[2]);
                    eles_num[0].nodes_id.push_back(nodes_num[3]);
                    
                    eles_num[1].nodes_id.push_back(nodes_num[1]);
                    eles_num[1].nodes_id.push_back(nodes_num[2]);
                    eles_num[1].nodes_id.push_back(nodes_num[3]);
                    eles_num[1].nodes_id.push_back(nodes_num[5]);
                    
                    eles_num[2].nodes_id.push_back(nodes_num[1]);
                    eles_num[2].nodes_id.push_back(nodes_num[3]);
                    eles_num[2].nodes_id.push_back(nodes_num[4]);
                    eles_num[2].nodes_id.push_back(nodes_num[5]);
                    //----------------------------------------------------------------
                    //Insert the number of elements to element vector
                    ele_temp.push_back(eles_num[0]);
                    ele_temp.push_back(eles_num[1]);
                    ele_temp.push_back(eles_num[2]);
                }
            }
        }
        
        nodes.push_back(nod_temp);
        eles.push_back(ele_temp);
    }
    
    return 1;
}

//---------------------------------------------------------------------------
//Transform angles into matrix
MathMatrix GenNetwork::Get_transformation_matrix(const double &sita, const double &pha)const
{
    //M = M_pha*M_sita
    //          |cos(pha) -sin(pha) 0|
    // M_pha  = |sin(pha)  cos(pha) 0|
    //          |   0         0     1|
    //
    //          | cos(sita)  0  sin(sita)|
    // M_sita = |     0      1      0    |
    //          |-sin(sita)  0  cos(sita)|
    //Calculate the matrix elements directly, instead of multiplying two matrices
    MathMatrix M(3,3);
    M.element[0][0] = cos(pha)*cos(sita);
    M.element[0][1] = -sin(pha);
    M.element[0][2] = cos(pha)*sin(sita);
    
    M.element[1][0] = sin(pha)*cos(sita);
    M.element[1][1] = cos(pha);
    M.element[1][2] = sin(pha)*sin(sita);
    
    M.element[2][0] = -sin(sita);
    M.element[2][2] = cos(sita);
    
    return M;
}

//---------------------------------------------------------------------------
//Calculate the angles of a verctor in the spherical coordinate
int GenNetwork::Get_angles_vector_in_spherial_coordinates(const Point_3D &normal, double &sita, double &pha)const
{
    if(normal.x==0&&normal.y==0&&normal.z==0) { hout << "Error, three elements of the vector are all zero!" << endl; return 0; }
    sita =  acos(normal.z/sqrt(normal.x*normal.x+normal.y*normal.y+normal.z*normal.z));
    if(normal.x==0&&normal.y==0) pha = 0;
    else if(normal.y>=0) pha = acos(normal.x/sqrt(normal.x*normal.x+normal.y*normal.y));
    else if(normal.y<0) pha = 2*PI - acos(normal.x/sqrt(normal.x*normal.x+normal.y*normal.y));
    
    return 1;
}

//---------------------------------------------------------------------------
//Calculate a group of equidistant points along the circumference which is on the plane defined by the center point of the circle and the normal vector
int GenNetwork::Get_points_circle_in_plane(const Point_3D &center, const double &trans_sita, const double &trans_pha, const double &radius, const int &num_sec, vector<Node> &nod_temp)const
{
    //Insert the center point firstly
    Node new_node(center.x, center.y, center.z);
    nod_temp.push_back(new_node);
    
    //Define the transformation matrix
    MathMatrix trans_mat(3,3);
    trans_mat = Get_transformation_matrix(trans_sita, trans_pha);
    
    //1D vector defined by a matrix
    MathMatrix Rvec(3,1);
    Rvec.element[0][0] = 0;
    Rvec.element[1][0] = 0;
    Rvec.element[2][0] = radius;
    
    //1D vector defined by a matrix
    MathMatrix Res(3,1);
    
    double sita, pha;
    sita = 0.5*PI;	//Defined on the XOY plane
    for(int i=0; i<num_sec; i++)
    {
        pha = i*2*PI/num_sec;
        MathMatrix matrix_temp = trans_mat*Get_transformation_matrix(sita, pha);
        Res = matrix_temp*Rvec;
        
        new_node.x = center.x + Res.element[0][0];
        new_node.y = center.y + Res.element[1][0];
        new_node.z = center.z + Res.element[2][0]; 
        
        //Insert the points on the circumference
        nod_temp.push_back(new_node);
    }
    
    return 1;
}

//---------------------------------------------------------------------------
//Calculate a group of projected points (which are on the plane with the center point of the circle and the normal vector) 
//which are projected from a group of points on the previous circumference and projected along the direction of line_vec
int GenNetwork::Get_projected_points_in_plane(const Point_3D &center, const Point_3D &normal, const Point_3D &line, const int &num_sec, vector<Node> &nod_temp)const
{
    //Record the total number of nodes after the previous generation
    const int nod_size = (int)nod_temp.size();  
    
    //Insert the center point
    Node new_node(center.x, center.y, center.z);
    nod_temp.push_back(new_node);
    
    const double vectors_dot_product = normal.x*line.x+normal.y*line.y+normal.z*line.z;
    
    if(vectors_dot_product==0.0) 
    {
        //Corresponding to three points: number 0, 1 and 2, the peak of this angle is at the point number 1. 
        hout << "Error: these two normal vectors are perpendicular to each other!" << endl;
        return 0; 
    }
    
    for(int i=num_sec; i>0; i--)
    {
        Point_3D point(center.x-nod_temp[nod_size-i].x, center.y-nod_temp[nod_size-i].y, center.z-nod_temp[nod_size-i].z);
        new_node.x = nod_temp[nod_size-i].x + (normal.x*point.x+normal.y*point.y+normal.z*point.z)*line.x/ vectors_dot_product;
        new_node.y = nod_temp[nod_size-i].y + (normal.x*point.x+normal.y*point.y+normal.z*point.z)*line.y/ vectors_dot_product;
        new_node.z = nod_temp[nod_size-i].z + (normal.x*point.x+normal.y*point.y+normal.z*point.z)*line.z/ vectors_dot_product;
        
        //Insert the points on the circumference
        nod_temp.push_back(new_node);
    }
    
    return 1;
}

//===========================================================================
