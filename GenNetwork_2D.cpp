//====================================================================================
//SOFTWARE:	Network of Eelectrically Conductive Nanocomposites (NECN)
//CODE FILE:	GenNetwork.cpp
//OBJECTIVE:	To generate networks with overlapping
//AUTHOR:		Fei Han; Angel Mora
//E-MAIL:			fei.han@kaust.edu.sa	;	angel.mora@kaust.edu.sa
//====================================================================================

#include "GenNetwork_2D.h"
#include "Geometry_2D.h"

//Generate 2D nanwire networks with ovelapping
int GenNetwork::Generate_nanowire_networks(const struct Geom_RVE &geom_rve, const struct Nanowire_Geo &nanowire_geo, vector<Point_2D> &cpoints, vector<double> &cnts_radius, vector<vector<int> > &cstructures)const
{
	// Define a vector of vectors for storing the cartesian coordinates of nanowires
    vector<vector<Point_2D> > cnts_points;
	//Use the Mersenne Twister for the random number generation
	if (Generate_network_threads_mt(geom_rve, nanowire_geo, cnts_points, cnts_radius)==0) return 0;
    //-----------------------------------------------------------------------------------------------------------------------------------------
    //Transform the 2D cnts_points into 1D cpoints and 2D cstructuers
    if (Transform_cnts_points(cnts_points, cpoints, cstructures)==0) return 0;	
	
	
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
    // Use random_device to generate a seed for Mersenne twister engine.
    std::random_device rd;
    // Use Mersenne twister engine to generate pseudo-random numbers.
    //Generate differnet engines for different variables
    std::mt19937 engine_x(rd());
    std::mt19937 engine_y(rd());
    std::mt19937 engine_pha(rd());
	
    // "Filter" MT's output to generate double values, uniformly distributed on the closed interval [0, 1].
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
	Rectangle gvcub(geom_rve.origin,geom_rve.len_x,geom_rve.wid_y);
	//generate a Rectangle to represent the extended domain
	Rectangle excub(geom_rve.ex_origin,geom_rve.ex_len,geom_rve.ey_wid);

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
        //Randomly generate the orientation of the nanowire 'Pha' is the orintation of the nanowire (-PI/2,PI/2)
        double cnt_pha;
        if(Get_uniform_direction_mt(nanowire_geo, cnt_pha, engine_pha, dist)==0) return 0;
		
        //---------------------------------------------------------------------------
        //The increased weight of each nanowire (If the different radii of nanowire are considered, the linear_density may be different in every nanowire)
        const double wei_para = nanowire_geo.linear_density;
		
        //---------------------------------------------------------------------------
        //Randomly generate a seed (initial point) of a CNT in the extended RVE	
        
        Point_2D cnt_poi;
        if(Get_seed_point_mt(excub, cnt_poi, engine_x, engine_y, dist)==0) return 0;
		
		int counter=1;
        //Check overlapping of the initial point 	
		 while(!Check_overlapping(Cnts_points, cnt_poi)) {
            if(Get_seed_point_mt(excub, cnt_poi, engine_x, engine_y, dist)==0) return 0;
            cnt_seed_count++;					//record the number of seed generations
            //hout << "Seed deleted" << endl;
            if (counter == MAX_ATTEMPTS) {
                hout << "Too many attempts to resolve overlapping of an intial CNT point (" << counter << " attempts). ";
                hout << cnt_poi.x << ' ' << cnt_poi.y << endl;
                return 0;
            }
            counter ++;
		} 
		
		
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
                if(Get_intersecting_point_RVE_edge(excub, new_cnt[0], cnt_poi, touch_edge)==0) {
                    hout << "Error in Generate_network_threads"<<endl;
                    return 0;
                }
                cnt_poi = touch_edge;
            }
            // At this point a nanowire is created with length 'cnt_length' and diameter '2*cnt_rad' and orientation 'cnt_pha'
	    // Calculate the accumulated area and the weight 
                double temp_length = Effective_length_given_region(gvcub, new_cnt[0],cnt_poi);
                if (temp_length > 0.0)
                {
                    area_sum += temp_length*2*cnt_rad;		//a accumulation on the area (Not the projection area)
                    wt_sum += temp_length*wei_para;	     	//a accumulation on the weight
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
    
/*     hout << "There were " << point_overlap_count_unique << " overlapping points and ";
    hout << point_overlap_count << " overlaps, " << endl; */
    
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
    if(nanowire_geo.dir_distrib_type=="random")
    {   
        //phai is chosen in [-PI/2, PI/2] with uniform distribution
        cnt_pha = -PI/2 + PI * dist(engine_pha);//*/
        
    }
	// ******************** This might have to be  changed ************
    else if(nanowire_geo.dir_distrib_type=="specific")
    {
		cnt_pha = -nanowire_geo.angle_max + 2*nanowire_geo.angle_max *dist(engine_pha);
    return 1;
}
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