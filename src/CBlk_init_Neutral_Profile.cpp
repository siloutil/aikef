
#include "CBlk.h"
#include "parameters.h"
#include "utils.h"
#include "absolute_Globals.h"

#include <iostream>
#include <fstream>
#include <math.h>
#include <cassert>

using namespace std;

extern D_REAL **delta_of_L;
extern D_REAL *CellVol_of_L;
extern D_REAL **Blk_Length_of;
extern INT32 i_j_k;


#if defined(nonadiabatic_gradPE_TERM) || defined(use_neutral_species_as_field)

//!-------------------------------------------------------------//
//! init_Neutral_Profile: 					//
//!-------------------------------------------------------------//

//! set neutral profile simply by using a new function
//! neutral_profile(NeutralSpecies, x[3], Nneutral, VNeutral[3])

void CBlock::which_neutral_profile(INT32 neutralSpecies, D_REAL *x, D_REAL *Nneutral, D_REAL *VNeutral)
{
	neutral_torus(neutralSpecies, x, Nneutral, VNeutral);
	// neutral_profile_Enceladus(neutralSpecies, x, Nneutral, VNeutral);
// 	init_Cometary_Neutral_Profile(neutralSpecies, x, Nneutral, VNeutral);
//   init_neutral_profile_test(neutralSpecies, x, Nneutral, VNeutral);
}	



void CBlock::init_Cometary_Neutral_Profile(INT32 neutralSpecies, D_REAL *x, D_REAL *Nneutral, D_REAL *VNeutral)
{
	
	//! Set Density
	if(vec_len2(x)>delta_of_L[MAX_LEVEL][0]*delta_of_L[MAX_LEVEL][0])
	{       
		Nneutral[0]     = COM_Q[neutralSpecies]*COM_RawDensity/vec_len2(x);
	}
	else
	{
		Nneutral[0]      = COM_Q[neutralSpecies]*COM_RawDensity/(delta_of_L[MAX_LEVEL][0]*delta_of_L[MAX_LEVEL][0]);
	}
	
	
	D_REAL abs_x = vec_len(x);
	if( abs_x < 1.e-2) 
		abs_x = 1.e-2;
				
	//! Set Radial Velocity
	VNeutral[0]       = x[0] * ((COM_v0/SI_v0) /abs_x);
	VNeutral[1]       = x[1] * ((COM_v0/SI_v0) /abs_x);
	VNeutral[2]       = x[2] * ((COM_v0/SI_v0) /abs_x);

}

//! ----------------------------------------------------------------------------------------
//! profile of neutral density
//! ----------------------------------------------------------------------------------------
void CBlock::neutral_profile_Enceladus(INT32 neutralSpecies, D_REAL *x,D_REAL *Nneutral,D_REAL *Vneutral){
	
	const D_REAL H_theta_gas = 1./((opening_angle_gas*M_PI/180.)*(opening_angle_gas*M_PI/180.));
	const D_REAL H_theta_dust = 1./((opening_angle_dust*M_PI/180.)*(opening_angle_dust*M_PI/180.));
	const D_REAL cosphi0 = cos(phi0);
	const D_REAL costheta0 = cos(theta0);
	const D_REAL costheta0_ChEx = cos(M_PI - theta0);
	const D_REAL sinphi0 = sin(phi0);
	const D_REAL sintheta0 = sin(theta0);
	
	D_REAL absr = vec_len(x);

	if(!use_test_scenario)
	{	
		if(absr*absr > R_Moon*R_Moon )
		{

			D_REAL nneutral=0;	    
			
	//			D_REAL axis[3] = {0,0,1};
			
			//! NOTE: BE CAREFUL ABOUT LENGTHS IN RE AND IN l0 !!!!!!!!!!!!!!!!!!!!!!!!!!!
	//			D_REAL pos_wrt_source[3] ={x_BlockNode[0]-RE*axis[0],x_BlockNode[1]-RE*axis[1],x_BlockNode[2]-RE*axis[2]};
			
	//			D_REAL betrag = vec_len(pos_wrt_source);
		
	//			D_REAL theta = acos( (axis[0]*pos_wrt_source[0] + axis[1]*pos_wrt_source[1] + axis[2]*pos_wrt_source[2])/(1.*betrag));
			
		//	nneutral = base_density/(1.*SI_n0)*RE*RE/(absr*absr)*exp(-(theta)*(theta)*H_theta_gas)*exp(-(absr-RE)/H_d);
			
		//          

			for(short source=0;source<8;source++)
			{
				
				
			
				//! positions of jets from spitale & porco in spherical coordinates
				double jets[9][2]={	  { 172*2*M_PI/360 ,  57*2*M_PI/360 },
							{ 169*2*M_PI/360 , 135*2*M_PI/360 },
							{ 171*2*M_PI/360 , 157*2*M_PI/360 },
							{ 163*2*M_PI/360 , 301*2*M_PI/360 },
							{ 170*2*M_PI/360 , 290*2*M_PI/360 }, //! alternative source IV from Dong 2011
							{ 169*2*M_PI/360 ,  18*2*M_PI/360 },
							{ 177*2*M_PI/360 , 219*2*M_PI/360 },
							{ 165*2*M_PI/360 ,  60*2*M_PI/360 },
							{ 172*2*M_PI/360 , 334*2*M_PI/360 } };


				D_REAL axis[3] = {cos(jets[source][1])*sin(jets[source][0]),sin(jets[source][1])*sin(jets[source][0]),cos(jets[source][0])};

					//! NOTE: BE CAREFUL ABOUT LENGTHS IN RE AND IN l0 !!!!!!!!!!!!!!!!!!!!!!!!!!!
				D_REAL pos_wrt_source[3] ={x[0]-R_Moon*axis[0],x[1]-R_Moon*axis[1],x[2]-R_Moon*axis[2]};
				
				D_REAL betrag = vec_len(pos_wrt_source);
			
				D_REAL theta = acos( (axis[0]*pos_wrt_source[0] + axis[1]*pos_wrt_source[1] + axis[2]*pos_wrt_source[2])/(1.*betrag));

			
	// 				//! neutral plume from Saur 2008 
	// 				//! ----------------------------------------------------------------------------------------------------------------
	// 	// 			nneutral = base_density/(1.*SI_n0)*RE*RE/(absr*absr)*exp(-(theta)*(theta)*H_theta_gas)*exp(-(absr-RE)/H_d);
	// 				//! ----------------------------------------------------------------------------------------------------------------
	// 
	// 				//! neutral plume from Dong 2011 
	// 				//! ----------------------------------------------------------------------------------------------------------------	
	// 			
				//! for E3
	// 				D_REAL Mach = 1.6;
	// 				D_REAL SourceStrength_SI[9] = {1.3e+22,1.3e+22,1.3e+22,0,0,1.3e+22,0,0} ;
	// 				//!NOTE: in 1/cm
	// 
				//! for E5
	// // 				D_REAL Mach = 1.4;
	// // 				D_REAL SourceStrength_SI[9] = {3.2e+22,3.2e+22,3.2e+22,0,0,0,3.2e+22,0,0} ;
	// 				//!NOTE: in 1/cm
	// 				
				//! for E7 
				D_REAL Mach = 2.2;
	// 				D_REAL SourceStrength_SI[9] = {1.e+22,1.e+22,1.e+22,3.5e+22,0.e+22,0,0,0,0} ;  //! Dong 2011
				D_REAL SourceStrength_SI[9] = {1.1e+22,1.1e+22,1.1e+22,0.e+22,1.e+22,0,0,0,0} ;//! NOTE: Best E7 from Dong2011
	//  				D_REAL Mach = 6.;
	//			D_REAL SourceStrength_SI[9] = {1.2e+22,1.2e+22,1.2e+22,0.7e+22,0,0,0,0,0} ;  //! Dong 2011
	//				D_REAL SourceStrength_SI[9] = {0.9e+22,0.9e+22,0.9e+22,0.2e+22,0.9e+22,0,0,0,0} ;//! NOTE: Best E7
	// 				//!NOTE: in 1/cm
	// 		 				 				 		
				nneutral += SourceStrength_SI[source]/(M_PI*betrag*betrag*SI_x0*SI_x0*1.e+4)*(
									2.*Mach*cos(theta)/(1.*sqrt(M_PI))*exp(-(Mach*Mach)) 
									+ exp(-(Mach*Mach*sin(theta)*sin(theta)))*(1.+2.*Mach*Mach*cos(theta)*cos(theta))
									*(1.+erf(Mach*cos(theta))));
	// 
	// 				//! ----------------------------------------------------------------------------------------------------------------
	// 								
			}
			//! NOTE: bei LOWRES fehlt am peak bissl was -> factor 2	
				

			D_REAL n0 = 3.e+05;
			D_REAL nbg = 3.3e+04;

			Nneutral[0] = 2.*(nneutral +  R_Moon*R_Moon/(absr*absr)*n0 + nbg)/(1.*SI_n0*1.e-06);  
			
			Vneutral[0] = 0;
			Vneutral[1] = 0;
			Vneutral[2] = 0;
		
		}
		else{
			Nneutral[0]=0;
			
			Vneutral[0] = 0;
			Vneutral[1] = 0;
			Vneutral[2] = 0;
		}
	}
	else//! test scenario
	{
		if(x[0]*x[0]+x[1]*x[1]<4*R_Moon*4*R_Moon && fabs(x[2])<2*R_Moon)
		Nneutral[0] = 1;
		else
		Nneutral[0] = 0;	
	}	
		
}



//!-------------------------------------------------------------//
//! for neutral density					//
//!-------------------------------------------------------------//
void CBlock::set_analytical_neutral_profile(INT32 neutralSpecies)
{
	INT32  ind[3];

	D_REAL *NeutralRho, *UnX, *UnY, *UnZ, *p, *beta_new;

	D_REAL x[3], Nneutral[1], Vneutral[3], cell_intern_r[3];
	
	memset(Nneutral,0,sizeof(D_REAL));
	memset(Vneutral,0,3*sizeof(D_REAL));
	memset(cell_intern_r,0,3*sizeof(D_REAL));

	//! Set pointer to intern (Block) field
	NeutralRho= Field_Type[id_numberdensity_neutralSpecies1 + neutralSpecies];
	UnX = Field_Type[id_velocity_neutralSpecies1 + neutralSpecies];
	UnY = UnX +num_nodes_in_block;
	UnZ = UnY +num_nodes_in_block;
	p  = Field_Type[id_pressure_neutralSpecies1 + neutralSpecies];
	beta_new = Field_Type[id_new_electron_beta_neutralSpecies1 + neutralSpecies];
	

	for (ind[0]=0; ind[0] < BlkNds_X; ind[0]++)
	 for (ind[1]=0; ind[1] < BlkNds_Y; ind[1]++)
	  for (ind[2]=0; ind[2] < BlkNds_Z; ind[2]++)
	  {
	
	
		i_j_k =  ind[0]*BlkNds_Y*BlkNds_Z
			+ind[1]*BlkNds_Z
			+ind[2];
		
			
		intern2normedCoords(x,cell_intern_r,ind);
		
		//! necessary if existing profile should be modified
		Nneutral[0] = NeutralRho[i_j_k];
		
		which_neutral_profile(neutralSpecies,x,Nneutral,Vneutral);
			
		NeutralRho[i_j_k] = Nneutral[0];	
		
		UnX[i_j_k] = Vneutral[0];
		UnY[i_j_k] = Vneutral[1];
		UnZ[i_j_k] = Vneutral[2];
		
		
		//!pressure is plasma beta at background density * normed density
		p[i_j_k]  = NeutralRho[i_j_k]*Neutral_Betas[neutralSpecies]; 
                                
		beta_new[i_j_k] = NewElectron_Betas[neutralSpecies];
	  }
	  
}	  

#endif



#if defined(use_ion_production_as_field)

//! decide which ion production profile is used	  
void CBlock::set_ion_production_profile(INT32 neutral_species)
{
	set_ion_production_profile_Enceladus(neutral_species);
}	


void CBlock::set_ion_production_profile_Enceladus(INT32 neutral_species)
{
	//! Saturn local time of all flybys
	D_REAL SLT[] = {22.6,17.0,17.1,23.2,22.7,22.5,22.5,10.9,3.7,9.3,3.9,3.6,8.8,8.8,23.8,23.8,23.8,0.6,0.6,0.5};
	//! Saturn's sub solar latitude of all flybys (in degrees)
	D_REAL SSL[] = {-22,-22,-21,-7.5,-5,-4.5,-4.25,1.275,1.565,3.954,4.132,5.543,7.132,7.435,11.39,11.62,12.05,13.8,13.92,14.28};

	D_REAL SLT_Degree, Sun_ENIS_angle, Sun_ENIS_vec[3];  
	  
	 SLT_Degree = SLT[flyby]/24.*360;
	 Sun_ENIS_angle = (270 - SLT_Degree)/360.*2.*M_PI;
	 
	 //! comment to SSL:
	 //! negative SSL is southern summer, that is a positive Sun_ENIS_vec[2]
	 //! thus pi + SSL with SSL<0 !
	 Sun_ENIS_vec[0] = cos(Sun_ENIS_angle)*sin(M_PI/2.+2.*M_PI*SSL[flyby]/360.);
	 Sun_ENIS_vec[1] = sin(Sun_ENIS_angle)*sin(M_PI/2.+2.*M_PI*SSL[flyby]/360.);
	 Sun_ENIS_vec[2] = cos(M_PI/2.+2.*M_PI*SSL[flyby]/360.);
	 
	bool Saturn_Shadow = false; 
	
	//! comment: length of geometric shadow in SLT:
	//! whole orbit is RS: 2 pi * 3.95
	//! divide by 24 to get 1 hour in RS or inverse to have 1 RS in hours
	//! ---> 1.934 local time hours
	//! no photionization from about 23 SLT to 1 SLT
	if(flyby==14 || flyby==15 || flyby==16)
	Saturn_Shadow = true;

	
	INT32  ind[3];

	D_REAL *IonProd, *UnX, *UnY, *UnZ;

	D_REAL x[3], Nneutral[1], Vneutral[3], cell_intern_r[3];
	
	memset(Nneutral,0,sizeof(D_REAL));
	memset(Vneutral,0,3*sizeof(D_REAL));
	memset(cell_intern_r,0,3*sizeof(D_REAL));

	//! Set pointer to intern (Block) field
	IonProd= Field_Type[id_density_ionProdSpecies1 + neutral_species];
	UnX = Field_Type[id_velocity_neutralSpecies1 + neutral_species];
	UnY = UnX +num_nodes_in_block;
	UnZ = UnY +num_nodes_in_block;

	
	if(!Saturn_Shadow)
	for (ind[0]=0; ind[0] < BlkNds_X; ind[0]++)
	 for (ind[1]=0; ind[1] < BlkNds_Y; ind[1]++)
	  for (ind[2]=0; ind[2] < BlkNds_Z; ind[2]++)
	  {
		i_j_k =  ind[0]*BlkNds_Y*BlkNds_Z
			+ind[1]*BlkNds_Z
			+ind[2];
		
		intern2normedCoords(x,cell_intern_r,ind);
			
		D_REAL temp[3];
		D_REAL x_scalar_Sun_ENIS_vec = x[0]*Sun_ENIS_vec[0]
						+ x[1]*Sun_ENIS_vec[1]
						+ x[2]*Sun_ENIS_vec[2];
		D_REAL abs2_Sun_ENIS_vec = vec_len2(Sun_ENIS_vec);
		
		temp[0] = x[0]-
			(Sun_ENIS_vec[0]*x_scalar_Sun_ENIS_vec)/abs2_Sun_ENIS_vec;
							
							
		temp[1] = x[1]-
			(Sun_ENIS_vec[1]*x_scalar_Sun_ENIS_vec)/abs2_Sun_ENIS_vec;

				
		temp[2] = x[2]-
			(Sun_ENIS_vec[2]*x_scalar_Sun_ENIS_vec)/abs2_Sun_ENIS_vec;
	
		
		D_REAL dist_to_shadow_axis = vec_len(temp);

		
		if(dist_to_shadow_axis < R_Obstacle && x_scalar_Sun_ENIS_vec >0)
		IonProd[i_j_k] = 0;

	  }
	  else
	  IonProd[i_j_k] = 0;
}

#endif

//!-------------------------------------------------------------//
//! set_RhoUi_extern: 								//
//!-------------------------------------------------------------//
void CBlock::set_RhoUi_extern(INT32 id_densityfield, INT32 id_velocityfield, short int* num_Nodes, FILE_REAL* Origin, FILE_REAL* Length, FILE_REAL* rho, INT32& num_values_not_in_extern_box)
{



	INT32  ind[3];
	INT32  num_extern_box_nodes;
	INT32 a,b,c, a_b_c;
	INT32 ap1_b_c, a_bp1_c, a_b_cp1;
	INT32 ap1_bp1_c, ap1_b_cp1, a_bp1_cp1, ap1_bp1_cp1;

	D_REAL *RHO, *UiX, *UiY, *UiZ;
	FILE_REAL *extern_UiX, *extern_UiY, *extern_UiZ;
	D_REAL r[3], extern_delta[3], shape_func[8];

	PARTICLE_REAL x_BlockNode[3], x_extern[3], cell_intern_r[3];
	memset(cell_intern_r,0,3*sizeof(PARTICLE_REAL));


	for(INT32 comp=0; comp<3; comp++)
	extern_delta[comp] = Length[comp]/(num_Nodes[comp]-1);


	//! Set pointer to extern field
	num_extern_box_nodes =   num_Nodes[0]
				*num_Nodes[1]
				*num_Nodes[2];

	extern_UiX = rho +num_extern_box_nodes;
	extern_UiY = extern_UiX +num_extern_box_nodes;
	extern_UiZ = extern_UiY +num_extern_box_nodes;
	
	
	//! Set pointer to intern (Block) field
	RHO = Field_Type[id_densityfield];

	UiX = Field_Type[id_velocityfield];
	UiY = UiX +num_nodes_in_block;
	UiZ = UiY +num_nodes_in_block;


	//! use a,b,c for extern mesh 
	//! use i,j,k for intern mesh
	for (ind[0]=0; ind[0] < BlkNds_X; ind[0]++)
	 for (ind[1]=0; ind[1] < BlkNds_Y; ind[1]++)
	  for (ind[2]=0; ind[2] < BlkNds_Z; ind[2]++)
	  {
	
	
		i_j_k =  ind[0]*BlkNds_Y*BlkNds_Z
			+ind[1]*BlkNds_Z
			+ind[2];
		
		//! get Coordinate of intern Block Node 
		intern2normedCoords(x_BlockNode, cell_intern_r, ind);
		
		//! find lower next neighbour node in read Mesh
		x_extern[0] = (x_BlockNode[0] +Origin[0])/extern_delta[0];
		x_extern[1] = (x_BlockNode[1] +Origin[1])/extern_delta[1];
		x_extern[2] = (x_BlockNode[2] +Origin[2])/extern_delta[2];

		 a  = int(x_extern[0]);
		 b  = int(x_extern[1]);
		 c  = int(x_extern[2]);


		if(   (x_extern[0]<0) || a>num_Nodes[0]-2
		    ||  x_extern[1]<0 || b>num_Nodes[1]-2
		    ||  x_extern[2]<0 || c>num_Nodes[2]-2)
		  {

			num_values_not_in_extern_box++;


		  }
		  else
		  {
		  
			r[0] = x_extern[0]-a;
			r[1] = x_extern[1]-b;
			r[2] = x_extern[2]-c;
			
			shape_func[0] = (1.-r[0])*(1.-r[1])*(1.-r[2]);
			shape_func[1] = (   r[0])*(1.-r[1])*(1.-r[2]);
			shape_func[2] = (1.-r[0])*(   r[1])*(1.-r[2]);
			shape_func[3] = (1.-r[0])*(1.-r[1])*(   r[2]);
			
			shape_func[4] = (   r[0])*(   r[1])*(1.-r[2]);
			shape_func[5] = (   r[0])*(1.-r[1])*(   r[2]);
			shape_func[6] = (1.-r[0])*(   r[1])*(   r[2]);
			shape_func[7] = (   r[0])*(   r[1])*(   r[2]);
			
			//! -----------------------------------------------
			a_b_c   =      a*num_Nodes[1]*num_Nodes[2]    +b*num_Nodes[2]      +c;
			
			//! -----------------------------------------------
			ap1_b_c = (a+1)*num_Nodes[1]*num_Nodes[2]     +b*num_Nodes[2]     +c;
			a_bp1_c =     a*num_Nodes[1]*num_Nodes[2] +(b+1)*num_Nodes[2]     +c;
			a_b_cp1 =     a*num_Nodes[1]*num_Nodes[2]     +b*num_Nodes[2] +(c+1);
			
			//! ------------------------------------------------
			ap1_bp1_c = (a+1)*num_Nodes[1]*num_Nodes[2] +(b+1)*num_Nodes[2]     +c;
			ap1_b_cp1 = (a+1)*num_Nodes[1]*num_Nodes[2]     +b*num_Nodes[2] +(c+1);
			a_bp1_cp1 =     a*num_Nodes[1]*num_Nodes[2] +(b+1)*num_Nodes[2] +(c+1);
			
			ap1_bp1_cp1 = (a+1)*num_Nodes[1]*num_Nodes[2] +(b+1)*num_Nodes[2] +(c+1);

			RHO[i_j_k]  =  rho[  a_b_c] * shape_func[0]
				      +rho[ap1_b_c] * shape_func[1]
				      +rho[a_bp1_c] * shape_func[2]
				      +rho[a_b_cp1] * shape_func[3]
			
				      +rho[  ap1_bp1_c] * shape_func[4]
				      +rho[  ap1_b_cp1] * shape_func[5]
				      +rho[  a_bp1_cp1] * shape_func[6]
				      +rho[ap1_bp1_cp1] * shape_func[7];


			UiX[i_j_k]  =  extern_UiX[  a_b_c] * shape_func[0]
				      +extern_UiX[ap1_b_c] * shape_func[1]
				      +extern_UiX[a_bp1_c] * shape_func[2]
				      +extern_UiX[a_b_cp1] * shape_func[3]
			
				      +extern_UiX[  ap1_bp1_c] * shape_func[4]
				      +extern_UiX[  ap1_b_cp1] * shape_func[5]
				      +extern_UiX[  a_bp1_cp1] * shape_func[6]
				      +extern_UiX[ap1_bp1_cp1] * shape_func[7];

			UiY[i_j_k]  =  extern_UiY[  a_b_c] * shape_func[0]
				      +extern_UiY[ap1_b_c] * shape_func[1]
				      +extern_UiY[a_bp1_c] * shape_func[2]
				      +extern_UiY[a_b_cp1] * shape_func[3]
			
				      +extern_UiY[  ap1_bp1_c] * shape_func[4]
				      +extern_UiY[  ap1_b_cp1] * shape_func[5]
				      +extern_UiY[  a_bp1_cp1] * shape_func[6]
				      +extern_UiY[ap1_bp1_cp1] * shape_func[7];

			UiZ[i_j_k]  =  extern_UiZ[  a_b_c] * shape_func[0]
				      +extern_UiZ[ap1_b_c] * shape_func[1]
				      +extern_UiZ[a_bp1_c] * shape_func[2]
				      +extern_UiZ[a_b_cp1] * shape_func[3]
			
				      +extern_UiZ[  ap1_bp1_c] * shape_func[4]
				      +extern_UiZ[  ap1_b_cp1] * shape_func[5]
				      +extern_UiZ[  a_bp1_cp1] * shape_func[6]
				      +extern_UiZ[ap1_bp1_cp1] * shape_func[7];
			
	      
 
			if(Flag[i_j_k])
			{
				RHO[i_j_k] = 0;
				UiX[i_j_k] = 0;
				UiY[i_j_k] = 0;
				UiZ[i_j_k] = 0;
			}	

		}
	  }

}


void CBlock::init_neutral_profile_test(INT32 neutralSpecies, D_REAL* x, D_REAL* Nneutral, D_REAL* VNeutral)
{
  //! Set density
  Nneutral[0]  = 2.00e8;
  
  //! Set velocity
  VNeutral[0] = V_sw[0];
  VNeutral[1] = V_sw[1];
  VNeutral[2] = V_sw[2];
}
// nrho[u_v_w]     = 2.00e9;
// ux[u_v_w]       = V_sw[0];
// uy[u_v_w]       = V_sw[1];
// uz[u_v_w]       = V_sw[2];
// p[u_v_w]        = 2.776e6;
// beta_new[u_v_w] = 1.38e-3;

void CBlock::neutral_torus(INT32 neutralSpecies, D_REAL* x, D_REAL* Nneutral, D_REAL* VNeutral)
{	
  	//! torus density with velocity zero
	VNeutral[0] = (VNeutral[0]*Nneutral[0])/(Nneutral[0]+neutral_torus_density);
	VNeutral[1] = (VNeutral[1]*Nneutral[0])/(Nneutral[0]+neutral_torus_density);
	VNeutral[1] = (VNeutral[2]*Nneutral[0])/(Nneutral[0]+neutral_torus_density);

	if(Nneutral[0]<neutral_torus_density)
	Nneutral[0] = neutral_torus_density;
}



