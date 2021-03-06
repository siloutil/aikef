

#include "CBlk.h"
#include "utils.h"
#include <math.h>
#include "parameters.h"
#include "absolute_Globals.h"
#include "gsl/gsl_rng.h"
#include "gsl/gsl_randist.h"


using namespace std;
extern D_REAL *dt_particle_of_L;
extern D_REAL **delta_of_L;

//!-----------------------------------------------------------------
//! D_REAL CBlock::checkChemicalReaction(
//!------------- ----------------------------------------------------
bool CBlock::checkChemicalReaction(D_REAL totalProcessRate,D_REAL dt_particle)
{
    //! beta = survival probability against charge exchange 
    D_REAL beta=totalProcessRate*dt_particle;
    
    //! Save Max Beta
    if(check_max_reaction_probability)
    if(beta > max_ChemicalReaction_Probability)
        max_ChemicalReaction_Probability = beta;
    
    if(beta > 1.*random()/RAND_MAX)
        return true;
    else
        return false;
}

//!-----------------------------------------------------------------
//! INT32 CBlock::determineDestinationSpecies(D_REAL totalProcessRate, D_REAL* processRates)
//!------------- ----------------------------------------------------
INT32 CBlock::determineDestinationSpecies(D_REAL totalProcessRate, D_REAL* processRates)
{
    //! First: determine intervalls according to fraction of reaction probability
    D_REAL dest_interval[num_Particle_Species+1];
    memset(dest_interval, 0,(num_Particle_Species+1)*sizeof(D_REAL));
    
    //! Intervall: 
    //!               0|-|---------|------|1
    //!dest_species     1   2         3
    for(INT32 dest_species=0; dest_species<num_Particle_Species+1; dest_species++) 
    {
        if(dest_species==0)
            dest_interval[dest_species] = processRates[dest_species]/totalProcessRate;
        else
            dest_interval[dest_species] = dest_interval[dest_species-1] + processRates[dest_species]/totalProcessRate;
    }
    
    //! Second: compare random number with reaction intervalls and decide
    //!             destination species (i.e which reaction takes place)
    D_REAL zufall = 1.*random()/RAND_MAX;
    INT32 counter=0;
    for(INT32 dest_species=0; dest_species<num_Particle_Species+1; dest_species++)
    {                                       
        if(zufall < dest_interval[dest_species])
        {
            return dest_species;
        }       
    }
    
    return -1;
}

//!-----------------------------------------------------------------
//! INT32 
//!------------- ----------------------------------------------------
INT32 CBlock::determineNeutralSpecies(INT32 dest_species, D_REAL* processRates, D_REAL* processRates_eachNeutralSpec)
{
    //! First: determine intervalls according to fraction of reaction probability
    D_REAL dest_interval[NUM_NEUTRAL_SPECIES];
    memset(dest_interval, 0,(NUM_NEUTRAL_SPECIES)*sizeof(D_REAL));

    //! Interval: 
    //!               0|-|---------|------|1
    //! neutral_species     1   2         3
    for(INT32 neutral_species=0; neutral_species<NUM_NEUTRAL_SPECIES; neutral_species++) 
    {
        if(neutral_species==0)
            dest_interval[neutral_species] = processRates_eachNeutralSpec[dest_species*NUM_NEUTRAL_SPECIES+neutral_species]/processRates[dest_species];
        else
            dest_interval[neutral_species] = dest_interval[neutral_species-1] + processRates_eachNeutralSpec[dest_species*NUM_NEUTRAL_SPECIES+neutral_species]/processRates[dest_species];
    }
    
    //! Second: compare random number with reaction intervalls and decide
    //!             destination species (i.e which reaction takes place)
    D_REAL zufall = 1.*random()/RAND_MAX;
    INT32 counter=0;
    for(INT32 neutral_species=0; neutral_species<NUM_NEUTRAL_SPECIES; neutral_species++) 
    {                                       
        if(zufall < dest_interval[neutral_species])
        {
            return neutral_species;
        }       
    }
    
    return -1;
}

//!-----------------------------------------------------------------
//! D_REAL CBlock::pickScalarFieldValue(INT32 field_id,particle* active_particle ,INT32 *ind)
//! Get value of a scalar field the particle position
//!------------- ----------------------------------------------------
D_REAL CBlock::pickScalarFieldValue(INT32 field_id,PARTICLE_REAL* vec_r ,INT32 *ind)
{
    D_REAL FieldValue = 0;
    D_REAL *Field_Cell;
    D_REAL shape_func[8];
        
    Field_Cell = Field_Type[field_id];
     
    //! indices for cell corners
    INT32 i_j_k   =  ind[0]*BlkNds_Y*BlkNds_Z    +  ind[1]*BlkNds_Z    +  ind[2];
    INT32 ip1_j_k = (ind[0]+1)*BlkNds_Y*BlkNds_Z     +ind[1]*BlkNds_Z     +ind[2];
    INT32 i_jp1_k =     ind[0]*BlkNds_Y*BlkNds_Z +(ind[1]+1)*BlkNds_Z     +ind[2];
    INT32 i_j_kp1 =     ind[0]*BlkNds_Y*BlkNds_Z     +ind[1]*BlkNds_Z +(ind[2]+1);
    INT32 ip1_jp1_k = (ind[0]+1)*BlkNds_Y*BlkNds_Z +(ind[1]+1)*BlkNds_Z     +ind[2];
    INT32 ip1_j_kp1 = (ind[0]+1)*BlkNds_Y*BlkNds_Z     +ind[1]*BlkNds_Z +(ind[2]+1);
    INT32 i_jp1_kp1 =     ind[0]*BlkNds_Y*BlkNds_Z +(ind[1]+1)*BlkNds_Z +(ind[2]+1);
    INT32 ip1_jp1_kp1 = (ind[0]+1)*BlkNds_Y*BlkNds_Z +(ind[1]+1)*BlkNds_Z +(ind[2]+1);
    
    
    //! --- now get shape function for interpolation ---
    shape_func[0] = (1.-vec_r[0])*(1.-vec_r[1])*(1.-vec_r[2]);
    shape_func[1] = (   vec_r[0])*(1.-vec_r[1])*(1.-vec_r[2]);
    shape_func[2] = (1.-vec_r[0])*(   vec_r[1])*(1.-vec_r[2]);
    shape_func[3] = (1.-vec_r[0])*(1.-vec_r[1])*(   vec_r[2]);
    
    shape_func[4] = (   vec_r[0])*(   vec_r[1])*(1.-vec_r[2]);
    shape_func[5] = (   vec_r[0])*(1.-vec_r[1])*(   vec_r[2]);
    shape_func[6] = (1.-vec_r[0])*(   vec_r[1])*(   vec_r[2]);
    shape_func[7] = (   vec_r[0])*(   vec_r[1])*(  vec_r[2]);
    
    //! --- INTERPOLATION ----------------------------------------
    FieldValue = Field_Cell[i_j_k]        * shape_func[0]
    + Field_Cell[ip1_j_k]      * shape_func[1]
    + Field_Cell[i_jp1_k]      * shape_func[2]
    + Field_Cell[i_j_kp1]      * shape_func[3]
    + Field_Cell[ip1_jp1_k]    * shape_func[4]
    + Field_Cell[ip1_j_kp1]    * shape_func[5]
    + Field_Cell[i_jp1_kp1]    * shape_func[6]
    + Field_Cell[ip1_jp1_kp1]  * shape_func[7];
    
    return FieldValue;
}

void CBlock::pickVectorFieldValue(INT32 field_id,PARTICLE_REAL* vec_r ,INT32 *ind, PARTICLE_REAL * vec_v)
{
        D_REAL *FieldX_Cell,*FieldY_Cell,*FieldZ_Cell;
        D_REAL shape_func[8];
        FieldX_Cell = Field_Type[field_id];
        FieldY_Cell = FieldX_Cell +num_nodes_in_block;
        FieldZ_Cell = FieldY_Cell +num_nodes_in_block;

        //! indices for cell corners
        INT32 i_j_k   =  ind[0]*BlkNds_Y*BlkNds_Z    +  ind[1]*BlkNds_Z    +  ind[2];
        INT32 ip1_j_k = (ind[0]+1)*BlkNds_Y*BlkNds_Z     +ind[1]*BlkNds_Z     +ind[2];
        INT32 i_jp1_k =     ind[0]*BlkNds_Y*BlkNds_Z +(ind[1]+1)*BlkNds_Z     +ind[2];
        INT32 i_j_kp1 =     ind[0]*BlkNds_Y*BlkNds_Z     +ind[1]*BlkNds_Z +(ind[2]+1);
        INT32 ip1_jp1_k = (ind[0]+1)*BlkNds_Y*BlkNds_Z +(ind[1]+1)*BlkNds_Z     +ind[2];
        INT32 ip1_j_kp1 = (ind[0]+1)*BlkNds_Y*BlkNds_Z     +ind[1]*BlkNds_Z +(ind[2]+1);
        INT32 i_jp1_kp1 =     ind[0]*BlkNds_Y*BlkNds_Z +(ind[1]+1)*BlkNds_Z +(ind[2]+1);
        INT32 ip1_jp1_kp1 = (ind[0]+1)*BlkNds_Y*BlkNds_Z +(ind[1]+1)*BlkNds_Z +(ind[2]+1);

        //     r[0] = vec_r[0];
        //     r[1] = vec_r[1];
        //     r[2] = vec_r[2];

        //! --- now get shape function for interpolation ---
        shape_func[0] = (1.-vec_r[0])*(1.-vec_r[1])*(1.-vec_r[2]);
        shape_func[1] = (   vec_r[0])*(1.-vec_r[1])*(1.-vec_r[2]);
        shape_func[2] = (1.-vec_r[0])*(   vec_r[1])*(1.-vec_r[2]);
        shape_func[3] = (1.-vec_r[0])*(1.-vec_r[1])*(   vec_r[2]);

        shape_func[4] = (   vec_r[0])*(   vec_r[1])*(1.-vec_r[2]);
        shape_func[5] = (   vec_r[0])*(1.-vec_r[1])*(   vec_r[2]);
        shape_func[6] = (1.-vec_r[0])*(   vec_r[1])*(   vec_r[2]);
        shape_func[7] = (   vec_r[0])*(   vec_r[1])*(  vec_r[2]);

        //! --- INTERPOLATION ----------------------------------------
        vec_v[0] = FieldX_Cell[i_j_k]        * shape_func[0]
                + FieldX_Cell[ip1_j_k]      * shape_func[1]
                + FieldX_Cell[i_jp1_k]      * shape_func[2]
                + FieldX_Cell[i_j_kp1]      * shape_func[3]
                + FieldX_Cell[ip1_jp1_k]    * shape_func[4]
                + FieldX_Cell[ip1_j_kp1]    * shape_func[5]
                + FieldX_Cell[i_jp1_kp1]    * shape_func[6]
                + FieldX_Cell[ip1_jp1_kp1]  * shape_func[7];

        vec_v[1] = FieldY_Cell[i_j_k]        * shape_func[0]
                + FieldY_Cell[ip1_j_k]      * shape_func[1]
                + FieldY_Cell[i_jp1_k]      * shape_func[2]
                + FieldY_Cell[i_j_kp1]      * shape_func[3]
                + FieldY_Cell[ip1_jp1_k]    * shape_func[4]
                + FieldY_Cell[ip1_j_kp1]    * shape_func[5]
                + FieldY_Cell[i_jp1_kp1]    * shape_func[6]
                + FieldY_Cell[ip1_jp1_kp1]  * shape_func[7];
        
        vec_v[2] = FieldZ_Cell[i_j_k]        * shape_func[0]
                + FieldZ_Cell[ip1_j_k]      * shape_func[1]
                + FieldZ_Cell[i_jp1_k]      * shape_func[2]
                + FieldZ_Cell[i_j_kp1]      * shape_func[3]
                + FieldZ_Cell[ip1_jp1_k]    * shape_func[4]
                + FieldZ_Cell[ip1_j_kp1]    * shape_func[5]
                + FieldZ_Cell[i_jp1_kp1]    * shape_func[6]
                + FieldZ_Cell[ip1_jp1_kp1]  * shape_func[7];
}

//!-------------------------------------------------------------//
//! charge exchange:                 //
//!-------------------------------------------------------------//
void CBlock::chemical_Reactions(INT64* num_exchanged_particle)
{
    
    particle* active_particle;
    //! cell indices
    INT32 i_j_k;
    INT32 ind[3];
    
    //! perform chemical reactions for each species if required
    for(INT32 species=0; species<num_Particle_Species; species++)
    if(species_does_chemical_reactions[species])	    
    {
        for(ind[0] = 1; ind[0] < BlkNds_X-1; ind[0]++)
            for(ind[1] = 1; ind[1] < BlkNds_Y-1; ind[1]++)
                for(ind[2] = 1; ind[2] < BlkNds_Z-1; ind[2]++)
                {
                    //! Calculate CellIndices
                    i_j_k   =  ind[0]*BlkNds_Y*BlkNds_Z
                            +  ind[1]*BlkNds_Z
                            +  ind[2];
                    //! for each particle in cell
                    for( INT32 part_index=0; part_index<num_MPiC[species][i_j_k];part_index++)
                    {
                        //! set pointer to particle at position part_index
                        active_particle = pArray[species][i_j_k] +part_index;
			
                        //! check if particle has already performed chemical reaction
                        //! This is done to avoid the following problem: 
                        //! 1) particle species 0 perform charge exchange to species 1
                        //! 2) after species 0  particle list loop starts for species 1 list 
                        //! 3) new species 1 particle should also perform chemical reaction ! PROBLEM!!
                        //! Hence check if particle is masked as already (v= v+10000); 
                        if(active_particle->v[0] < 1000.)
                        {
                            //! Calculate Process Rates [1/s](SI)
                            D_REAL totalProcessRate = 0.;
                            D_REAL processRates[num_Particle_Species+1];
                            memset(processRates,0,(num_Particle_Species+1)*sizeof(D_REAL));
			    
		            D_REAL processRates_eachNeutralSpec[(num_Particle_Species+1)*(NUM_NEUTRAL_SPECIES)];
                            memset(processRates_eachNeutralSpec,0,(num_Particle_Species+1)*NUM_NEUTRAL_SPECIES*sizeof(D_REAL));
			    			    			            
			    D_REAL NeutralDensity[NUM_NEUTRAL_SPECIES];
				
			    if(NUM_NEUTRAL_SPECIES>0)
			    for(INT32 neutral_species=0;neutral_species<NUM_NEUTRAL_SPECIES;neutral_species++)
				NeutralDensity[neutral_species] = getNeutralDensityAtPoint(active_particle->rel_r,ind,neutral_species);	
			    
			    
                            //! Set process Rates
                            for(INT32 dest_species=0;dest_species<num_Particle_Species+1;dest_species++)
                            {
                                if(dest_species<num_Particle_Species)
                                {				    
					processRates[dest_species]=0;
					
					for(INT32 neutral_species=0; neutral_species<NUM_NEUTRAL_SPECIES; neutral_species++)
					{						
						D_REAL vel_dependent_rate=calc_reaction_rate(active_particle->v,species,dest_species,neutral_species);
						
						processRates_eachNeutralSpec[dest_species*NUM_NEUTRAL_SPECIES+neutral_species] = vel_dependent_rate*NeutralDensity[neutral_species];
						processRates[dest_species] += vel_dependent_rate*NeutralDensity[neutral_species];
					}						
				
					totalProcessRate += processRates[dest_species];
                                }
                                else
                                {
				    D_REAL RecombinationRate = 0 ;
				    D_REAL ElectronDensity = 0;
					
				    if(Recombination_for_Species[species])
				    {	    
					//! Recombination
					//! Reaction H3O+ + e -> H2O + H (Puhl-Quinn 1995)
					RecombinationRate = pickScalarFieldValue(id_recomb_Species1+species,active_particle->rel_r,ind);
					
					//! Total Electron Denstiy [1/cm^3]
					ElectronDensity = pickScalarFieldValue(id_rho_np1,active_particle->rel_r,ind);
					
				    }
				    
                                    //!Recombination
                                    processRates[dest_species] = RecombinationRate*ElectronDensity; 
                                    totalProcessRate += RecombinationRate*ElectronDensity;
  
                                }
                            } // end set process Rates
                            
                            //! Check for max Rates
                            if(check_max_reaction_rates)
                            for(INT32 dest_species=0; dest_species<num_Particle_Species+1;dest_species++)
                            {
                                if(processRates[dest_species]>max_ChemicalReaction_Rates[species*(num_Particle_Species+1)+dest_species])
                                    max_ChemicalReaction_Rates[species*(num_Particle_Species+1)+dest_species] = processRates[dest_species];
                            }
                            
                            //! Check if Process takes place
                            if(checkChemicalReaction(totalProcessRate,dt_particle_of_L[RLevel]))
                            {
                                //! A Chemical Reaction takes place
                                
                                //! if chargeexchange occurs, mark particle fortracing
                                #ifdef TRACK_PARTICLE
                                active_particle->exchange++;
                                #endif
                                
                               
                                //! calculate reaction process
                                INT32 destinationSpecies = determineDestinationSpecies(totalProcessRate, processRates);
                                if(destinationSpecies==-1)
                                {
                                    log_file << "ERROR: destination Species wrong "<<endl;
                                    destinationSpecies = species;
                                }   
                                
                                if(destinationSpecies<num_Particle_Species)
				{	
					//! in principle, it is possible that the identical source and destination ions occur for reactions with
					//! different neutral species, e.g. the hypothetical reactions
					//! H+ + H2O -> OH+ + H2
					//! H+ + CO  -> OH+ + C
					//! Therefore, neutral species is determined in the same manner as the destination species
					INT32 NeutralSpecies = determineNeutralSpecies(destinationSpecies,processRates,processRates_eachNeutralSpec);    
					
					//! Set new particle velocity according to neutral species of reaction

					//! first neutral bulk velocity
					getNeutralVelocityAtPoint(active_particle->rel_r,ind, NeutralSpecies,active_particle->v);
					
					//! second thermal velocity
					PARTICLE_REAL vth[3];
					//! it is important to initialize vth, oth. error in case 
					//! temperature of species is zero.
					memset(vth, 0 ,3*sizeof(PARTICLE_REAL));
					
					vth[0]= gsl_ran_gaussian_ziggurat(randGen_of_species[species], Neutral_vth[NeutralSpecies]);
					vth[1]= gsl_ran_gaussian_ziggurat(randGen_of_species[species], Neutral_vth[NeutralSpecies]);
					vth[2]= gsl_ran_gaussian_ziggurat(randGen_of_species[species], Neutral_vth[NeutralSpecies]);	
					
					active_particle->v[0] += vth[0];
					active_particle->v[1] += vth[1];
					active_particle->v[2] += vth[2];
					
				}
				//! mark particle by adding a unusual high velocity. This is neccessary, 
				//! since otherwise one particle could undergo more than one reaction
				//! in one timestep (i. e. change species more often and therefore
				//! cause errors in chemistry)
				//! this has to be done before particle changes arrays
				//! BUT AFTER NEW VELOCITY IS SET
				active_particle->v[0] += 10000.;
								
                                //! change particle's list if neccessary
                                if(destinationSpecies != species)
                                {
                                    //! Check If ChargeExchange and IonNeutral Collision or Recombination
                                    if(destinationSpecies<num_Particle_Species)
                                    {

					
                                        //! ChargeExchange and IonNeutral Collision
                                        //! add particle to list of last (obstacle)ion species
                                        //! num_MPiC[destinationSpecies][i_j_k] is increased in add_particle_to_pArray function
                                        add_particle_to_pArray(destinationSpecies, i_j_k,active_particle);
        
                                    }
                                    else
				    {
				      //! Recombination 
				      
				      //! update statistics
				      num_total_particles--;
				      num_total_particles_in_L[RLevel]--;
				    }
				    
                                    //! DELETE ACTIVE PARTICLE from list of species
                                    //! move next up to last particle of list on top active_particle
                                    memmove(pArray[species][i_j_k] +(part_index),
                                            pArray[species][i_j_k] +(part_index+1),
                                            (num_MPiC[species][i_j_k] -(part_index+1))*sizeof(particle));
                                    
                                    //! update statistics
                                    num_MPiC[species][i_j_k]--;
                                    
                                    //! part_index=-1 may temporary occcur here
                                    part_index--;
                                    
                                }// end destinationSpecies != species
                                
                                //! for statistics
                                num_exchanged_particle[species*(num_Particle_Species+1)+destinationSpecies]++;
                                
                            } // end check Chemical Reaction
                        }// end check particle masked (v=1000)
                    } // end for each particle in cell
                } // end for BlkNds Z,Y,Z
    } // end for species
    
    //! unMark particles as not charge-exchanged (v>1000)
    for(INT32 species=0; species<num_Particle_Species; species++)
    {
        
        for(ind[0] = 1; ind[0] < BlkNds_X-1; ind[0]++)
            for(ind[1] = 1; ind[1] < BlkNds_Y-1; ind[1]++)
                for(ind[2] = 1; ind[2] < BlkNds_Z-1; ind[2]++)
                {
                    
                    i_j_k   =  ind[0]*BlkNds_Y*BlkNds_Z
                    +ind[1]*BlkNds_Z
                    +ind[2];
                    
                    
                    for(INT32 part_index=0; part_index<num_MPiC[species][i_j_k]; part_index++)
                    {
                        active_particle = pArray[species][i_j_k] +part_index;
                        
                        if(active_particle->v[0] > 1000.)
                            active_particle->v[0] += -10000.;
                    }
                }
    }//! end demark particles
}

//! 
//! prepare_Recombination_Density
//!

//! Function copy Field from total density to recombined density

void CBlock::prepare_Recombination_Density()
{
        //! Delete first old values
        set_zero_Field(id_rho_np1_recombined);
        //! Copy total density to recombined density field
        copy_Field(id_rho_np1_recombined,id_rho_np1);
}

//! 
//! calc_Recombination_Density
//!

//! function calculate the recombined density 

void CBlock::calc_Recombination_Density()
{
        //! substract density gathered after from density gathered before recombination process 
        add_multipliedField(id_rho_np1_recombined, id_rho_np1, -1.);
}


//!-------------------------------------------------------------//
//! calc_RecombinationAlphaField
//! Use Eq. 11 after Gombosi 1996
//!-------------------------------------------------------------//
void CBlock::calc_RecombinationAlphaField(void)
{	
	D_REAL *RecombinationRate;
	D_REAL *Te;
    
	Te = Field_Type[id_ElectronTemperature];
	
	for(INT32 species=0; species<num_Particle_Species; species++)
	if(Recombination_for_Species[species])	
	{	
		RecombinationRate = Field_Type[id_recomb_Species1 +species];
		
		for (INT32 node=0; node < num_nodes_in_block; node++)
		{
			if(Te[node]<=200)
			{
			RecombinationRate[node] = sqrt(300/Te[node]);
			}
			else 
			{
			D_REAL exp=0.2553-0.1633*log(Te[node]);
			RecombinationRate[node] = 2.342*pow(Te[node],exp);
			}
		}
	}	
}

//!-----------------------------------------------------------------
//! calculate the reaction rate for all reaction of species
//! if already specified in parameter.cpp, just set values
//! else calculate rate for all reactions of species
//! according to particle's velocity
//!-----------------------------------------------------------------
D_REAL calc_reaction_rate(PARTICLE_REAL* v, INT32 species, INT32 dest_species, INT32 neutral_species)
{
	
	D_REAL rate_of_reaction=ReactionRate[species][dest_species][neutral_species];
	
	if(use_velocity_dependent_reaction_rates==false)
		return rate_of_reaction;
	else
	{
		//! to speed up calculation of velocity dependent rates, these rates are calculated once at beginning of simulation
		//! and binned into array velocity_rate
		
		INT32 v_array = int(vec_len(v)*SI_v0*1e-03);
		if(v_array<1)
		v_array=1;
		if(v_array>=30)
		v_array=29;
		
		//! rates for charge_exchange channel of O+ + H2O -> H2O+ + O	
		if(species==0)
		{
			if(dest_species == 2)
			rate_of_reaction += velocity_rate[v_array][species];
		}
		
		//! rates for OH+ + H2O -> H2O+ + OH	
		//! (identical to O+ + H2O, just with different mass)
		if(species==1)
		{
			if(dest_species == 2)
			rate_of_reaction += velocity_rate[v_array][species];
		}
		
		//! rates for H2O+ + H2O -> H2O+ + H2O	
		if(species==2)
		{
			if(dest_species == 2)
			rate_of_reaction += velocity_rate[v_array][species];
		}
		
		//! rates for charge exchange channel H+ + H2O -> H2O+ + H	
		if(species==4)
		{
			if(dest_species == 2)
			rate_of_reaction += velocity_rate[v_array][species];
		}
		
		return rate_of_reaction;

	}

}


