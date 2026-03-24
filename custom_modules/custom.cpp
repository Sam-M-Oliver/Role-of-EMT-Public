/*
###############################################################################
# If you use PhysiCell in your project, please cite PhysiCell and the version #
# number, such as below:                                                      #
#                                                                             #
# We implemented and solved the model using PhysiCell (Version x.y.z) [1].    #
#                                                                             #
# [1] A Ghaffarizadeh, R Heiland, SH Friedman, SM Mumenthaler, and P Macklin, #
#     PhysiCell: an Open Source Physics-Based Cell Simulator for Multicellu-  #
#     lar Systems, PLoS Comput. Biol. 14(2): e1005991, 2018                   #
#     DOI: 10.1371/journal.pcbi.1005991                                       #
#                                                                             #
# See VERSION.txt or call get_PhysiCell_version() to get the current version  #
#     x.y.z. Call display_citations() to get detailed information on all cite-#
#     able software used in your PhysiCell application.                       #
#                                                                             #
# Because PhysiCell extensively uses BioFVM, we suggest you also cite BioFVM  #
#     as below:                                                               #
#                                                                             #
# We implemented and solved the model using PhysiCell (Version x.y.z) [1],    #
# with BioFVM [2] to solve the transport equations.                           #
#                                                                             #
# [1] A Ghaffarizadeh, R Heiland, SH Friedman, SM Mumenthaler, and P Macklin, #
#     PhysiCell: an Open Source Physics-Based Cell Simulator for Multicellu-  #
#     lar Systems, PLoS Comput. Biol. 14(2): e1005991, 2018                   #
#     DOI: 10.1371/journal.pcbi.1005991                                       #
#                                                                             #
# [2] A Ghaffarizadeh, SH Friedman, and P Macklin, BioFVM: an efficient para- #
#     llelized diffusive transport solver for 3-D biological simulations,     #
#     Bioinformatics 32(8): 1256-8, 2016. DOI: 10.1093/bioinformatics/btv730  #
#                                                                             #
###############################################################################
#                                                                             #
# BSD 3-Clause License (see https://opensource.org/licenses/BSD-3-Clause)     #
#                                                                             #
# Copyright (c) 2015-2021, Paul Macklin and the PhysiCell Project             #
# All rights reserved.                                                        #
#                                                                             #
# Redistribution and use in source and binary forms, with or without          #
# modification, are permitted provided that the following conditions are met: #
#                                                                             #
# 1. Redistributions of source code must retain the above copyright notice,   #
# this list of conditions and the following disclaimer.                       #
#                                                                             #
# 2. Redistributions in binary form must reproduce the above copyright        #
# notice, this list of conditions and the following disclaimer in the         #
# documentation and/or other materials provided with the distribution.        #
#                                                                             #
# 3. Neither the name of the copyright holder nor the names of its            #
# contributors may be used to endorse or promote products derived from this   #
# software without specific prior written permission.                         #
#                                                                             #
# THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" #
# AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE   #
# IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE  #
# ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE   #
# LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR         #
# CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF        #
# SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS    #
# INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN     #
# CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)     #
# ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE  #
# POSSIBILITY OF SUCH DAMAGE.                                                 #
#                                                                             #
###############################################################################
*/

#include "./custom.h"

void create_cell_types( void )
{
	// set the random seed 
	SeedRandom( parameters.ints("random_seed") );  
	
	/* 
	   Put any modifications to default cell definition here if you 
	   want to have "inherited" by other cell types. 
	   
	   This is a good place to set default functions. 
	*/ 
	
	initialize_default_cell_definition(); 
	cell_defaults.phenotype.secretion.sync_to_microenvironment( &microenvironment ); 
	
	cell_defaults.functions.volume_update_function = standard_volume_update_function;
	cell_defaults.functions.update_velocity = standard_update_cell_velocity;

	cell_defaults.functions.update_migration_bias = NULL; 
	cell_defaults.functions.update_phenotype = NULL; // update_cell_and_death_parameters_O2_based; 
	cell_defaults.functions.custom_cell_rule = NULL; 
	cell_defaults.functions.contact_function = NULL; 
	
	cell_defaults.functions.add_cell_basement_membrane_interactions = NULL; 
	cell_defaults.functions.calculate_distance_to_membrane = NULL; 
	
	/*
	   This parses the cell definitions in the XML config file. 
	*/
	
	initialize_cell_definitions_from_pugixml(); 

	/*
	   This builds the map of cell definitions and summarizes the setup. 
	*/
		
	build_cell_definitions_maps(); 

	/*
	   This intializes cell signal and response dictionaries 
	*/

	setup_signal_behavior_dictionaries(); 	

	/*
       Cell rule definitions 
	*/

	setup_cell_rules(); 

	/* 
	   Put any modifications to individual cell definitions here. 
	   
	   This is a good place to set custom functions. 
	*/ 
	
	cell_defaults.functions.update_phenotype = phenotype_function; 
	cell_defaults.functions.custom_cell_rule = custom_function; 
	cell_defaults.functions.contact_function = contact_function; 
	
	Cell_Definition *pC = find_cell_definition("cancer");
	pC->functions.update_phenotype = cancer_phenotype;
	
	/*
	   This builds the map of cell definitions and summarizes the setup. 
	*/
		
	display_cell_definitions( std::cout ); 
	
	return; 
}

void setup_microenvironment( void )
{
	// set domain parameters 
	
	// put any custom code to set non-homogeneous initial conditions or 
	// extra Dirichlet nodes here. 
	
	// initialize BioFVM 
	
	initialize_microenvironment(); 	
	
	return; 
}

void setup_tissue( void )
{
	double Xmin = microenvironment.mesh.bounding_box[0]; 
	double Ymin = microenvironment.mesh.bounding_box[1]; 
	double Zmin = microenvironment.mesh.bounding_box[2]; 

	double Xmax = microenvironment.mesh.bounding_box[3]; 
	double Ymax = microenvironment.mesh.bounding_box[4]; 
	double Zmax = microenvironment.mesh.bounding_box[5]; 
	
	if( default_microenvironment_options.simulate_2D == true )
	{
		Zmin = 0.0; 
		Zmax = 0.0; 
	}
	
	double Xrange = Xmax - Xmin; 
	double Yrange = Ymax - Ymin; 
	double Zrange = Zmax - Zmin; 
	
	// create some of each type of cell 
	
	Cell* pC;
	
	double max_distance = parameters.doubles("max_initial_distance");
	
	Cell_Definition* pCD1 = find_cell_definition( "cancer" );
	
	for( int k=0 ; k < parameters.ints( "number_of_cells" ); k++ )
	{
		std::vector<double> position = {0,0,0};
		double r = cbrt(UniformRandom())* max_distance;
		double theta = UniformRandom() * 6.2831853;
		double rho = UniformRandom() * 6.2831853;
		position[0] = r*cos(theta)*sin(rho);
		position[1] = r*sin(theta)*sin(rho);
		position[2] = r*cos(rho);
		pC = create_cell( *pCD1 );
		pC->assign_position( position );
	}
	
	//placing ovcar cells
	if (parameters.ints("cancer_type")==0)
	{
		std::cout << "Placing cells of type " << "Ovcar" << " ... " << std::endl;
	}
	
	//placing skov cells
	if (parameters.ints("cancer_type")==1)
	{
		std::cout << "Placing cells of type " << "Skov" << " ... " << std::endl;
	}
	std::cout << std::endl;
	
	// load cells from your CSV file (if enabled)
	load_cells_from_pugixml(); 	
	
	return; 
}

std::vector<std::string> my_coloring_function( Cell* pCell )
{ return paint_by_number_cell_coloring(pCell); }

void phenotype_function( Cell* pCell, Phenotype& phenotype, double dt )
{ return; }

void custom_function( Cell* pCell, Phenotype& phenotype , double dt )
{ return; } 

void contact_function( Cell* pMe, Phenotype& phenoMe , Cell* pOther, Phenotype& phenoOther , double dt )
{ return; } 

void cancer_phenotype( Cell *pC, Phenotype &p, double dt )
{
	Cell_Definition *pCD = find_cell_definition("cancer");
		
	//finding substrate indicies
	static int nO2 = microenvironment.find_density_index("oxygen");
	static int npac = microenvironment.find_density_index("drug");
	static int nsig = microenvironment.find_density_index("signal");
	
	//finding phenotype indicies
	static int nNecro = p.death.find_death_model_index("Necrosis");
	static int nApop = p.death.find_death_model_index("Apoptosis");
	static int nC = pC->custom_data.find_variable_index("cadherin_level"); 
	
	//finding parameter values
	double pO2 = pC->nearest_density_vector()[nO2];
	double psig = pC->nearest_density_vector()[nsig];
	double ppac = pC->nearest_density_vector()[npac];
	double ppres = pC->state.simple_pressure; 
	double base_apop_rate = pCD->phenotype.death.rates[nApop];
	double base_prof_rate = pCD->phenotype.cycle.data.transition_rate(0, 1);
	
	//variables here mainly used for future sensitivity analysis
	double base_prob = 0.0;
	double apop_rate = 0.0;
	double cycle_rate = 0.0;
	
	//Hill functions for 6 main variables in each cell (cycling, emt rates)
	double emt_oxygen_rate = pC->custom_data["pO2_transition_impact"]/10.0 - pC->custom_data["pO2_transition_impact"]/10.0 * pow(pO2, 4)/(pow((pC->custom_data["pO2_transition_saturation"]), 4)+pow(pO2, 4));
	double emt_signal_on = pC->custom_data["psig_transition_impact"] * pow(psig, 4)/(pow((pC->custom_data["bystander_threshold"]), 4)+pow(psig, 4));
	double emt_base_rate = pC->custom_data["base_ovcar_prob"] * pow(pC->custom_data[nC], 4)/(pow(6.5, 4)+pow(pC->custom_data[nC], 4));
	double cycling_pressure = pC->custom_data["ppres_cycling_impact"] - pC->custom_data["ppres_cycling_impact"] * pow(ppres, 4)/(pow((pC->custom_data["ppres_cycling_threshold"]), 4)+pow(ppres, 4));
	double cycling_O2 = pC->custom_data["pO2_cycling_impact"] * pow(pO2, 4)/(pow(pC->custom_data["pO2_cycling_saturation"], 4)+pow(pO2, 4));
	double cycling_base_rate = pC->custom_data["base_ovcar_rate"] - pC->custom_data["base_ovcar_rate"] * pow(pC->custom_data[nC], 4)/(pow(6.5, 4)+pow(pC->custom_data[nC], 4));

	
	if (p.death.dead == true)
	{
		p.secretion.set_all_secretion_to_zero();
		p.secretion.set_all_uptake_to_zero();
		pC->functions.update_phenotype = NULL;
	}
	
	// ovcar has standard death rates and cycle rates
	if (parameters.ints("cancer_type") ==0) //ovcar
	{	
		apop_rate = base_apop_rate;
		cycle_rate = base_prof_rate * ( 1 + cycling_pressure + cycling_O2 + cycling_base_rate ) * fmax(0,1-((*all_cells).size()/parameters.doubles("carrying_capacity")));
	}
		
	// skov has higher emt rates, lower cycling and death rates
	if (parameters.ints("cancer_type") ==1) //skov
	{
		emt_oxygen_rate = 5 * emt_oxygen_rate;
		emt_signal_on = 5 * emt_signal_on;
		emt_base_rate = 5 * emt_base_rate;
		apop_rate = base_apop_rate/4.0;
		cycle_rate = base_prof_rate * ( 1 + cycling_pressure + cycling_O2 + cycling_base_rate ) * fmax(0,1-((*all_cells).size()/parameters.doubles("carrying_capacity"))) / 4.0;
			
		//if skov and on the edge of the tumour, cells are epithelial
		if (pO2 > pC->custom_data["pO2_transition_threshold"])
		{
			pC->custom_data[nC] = 0;
		}
	}
	
	double jump_up_prob = (base_prob+emt_oxygen_rate+emt_signal_on-base_prob*emt_oxygen_rate-base_prob*emt_signal_on-emt_oxygen_rate*emt_signal_on+base_prob*emt_oxygen_rate*emt_signal_on);
	double jump_down_prob = parameters.doubles("met_on")*((pC->custom_data["base_ovcar_prob"]-base_prob)+(pC->custom_data["pO2_transition_impact"]-emt_oxygen_rate)-(pC->custom_data["base_ovcar_prob"]-base_prob)*(pC->custom_data["pO2_transition_impact"]-emt_oxygen_rate));
			
	//created these jump up (or jump down probabilities if met is included)
	if (pC->custom_data[nC] == 0 )
	{
		double r=UniformRandom();
		if (r < jump_up_prob)
		{
			pC->custom_data[nC] = pC->custom_data[nC] + 1;
		}
	}

	else if (pC->custom_data[nC] > 0.5 && pC->custom_data[nC] < 12.5 )
	{
		double r=UniformRandom();
		if (r < jump_up_prob)
		{
			pC->custom_data[nC] = pC->custom_data[nC] + 1;
		}
		else if (r < jump_up_prob + jump_down_prob)
		{
			pC->custom_data[nC] = pC->custom_data[nC] - 1;
		}
	}
	
	else if (pC->custom_data[nC] == 13 )
	{
		double r=UniformRandom();
		if (r < jump_down_prob)
		{
			pC->custom_data[nC] = pC->custom_data[nC] - 1;
		}
	}
		
		//update the phentotypic behaviours
		p.secretion.secretion_rates[nsig] = 50*pow(pC->custom_data[nC], 4)/(pow(6.5, 4)+pow(pC->custom_data[nC], 4));
		
		p.death.rates[nApop] = apop_rate;
		p.cycle.data.transition_rate(0, 1) = cycle_rate * fmax(0,1-((*all_cells).size()/parameters.doubles("carrying_capacity")));
}

//impliment the drup
void drug_administration()
	{
		bool drug_on = false;
		static int npac = microenvironment.find_density_index("drug");
		
		// no dirichlet conditions if not at treatment time
		
		// dirichlet conditions if at treatment time
		
		if ((PhysiCell_globals.current_time > parameters.doubles("drug_start_time")*60) && (((parameters.doubles("drug_start_time")+1)*60 ) > PhysiCell_globals.current_time ))
		{
			drug_on = true;
			microenvironment.set_substrate_dirichlet_activation( npac , drug_on );
		}
		
		else 
		{
			drug_on = false;
			microenvironment.set_substrate_dirichlet_activation( npac , drug_on );
		}
	}

std::vector<std::string> custom_coloring_function( Cell* pC )
{
	
	std::vector<std::string> output = paint_by_number_cell_coloring(pC);
	// get cell definitions
	static Cell_Definition* pCD = find_cell_definition( "cancer" );
	int nC = pC->custom_data.find_variable_index( "cadherin_level" );
	// get cadherin rating
	double s = pC->custom_data[nC];
	// make color based on cadherin rating
	int color = (int) round( 19.6 * s );
	char szColor [1024];
	
	sprintf( szColor, "rgb(%u,%u,%u)", color, 255-color, 0);
	// modify output
	output[0] = szColor;
	output[2] = szColor;
	output[3] = szColor;
	return output;
		
}