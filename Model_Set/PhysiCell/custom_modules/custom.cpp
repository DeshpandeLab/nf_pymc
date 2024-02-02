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

extern std::vector<int> InitialFraction;

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
	   Put any modifications to individual cell definitions here.

	   This is a good place to set custom functions.
	*/

	cell_defaults.functions.update_phenotype = phenotype_function;
	cell_defaults.functions.custom_cell_rule = custom_function;
	cell_defaults.functions.contact_function = contact_function;

	Cell_Definition* pCD = find_cell_definition( "cancer");
	pCD->functions.update_phenotype = cancer_phenotype;
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

	Cell* pCell = NULL;
	Cell_Definition* pCD = find_cell_definition( "cancer" );
	double cell_radius = cell_defaults.phenotype.geometry.radius;
	double cell_spacing = 0.95 * 2.0 * cell_radius;
	double tumor_radius = parameters.doubles("tumor_radius");

	double x = 0.0;
	double x_outer = tumor_radius;
	double y = 0.0;

	if( default_microenvironment_options.simulate_2D == true ){
        int n = 0;
        while( y < tumor_radius )
        {
            x = 0.0;
            if( n % 2 == 1 )
            { x = 0.5*cell_spacing; }
            x_outer = sqrt( tumor_radius*tumor_radius - y*y );

            while( x < x_outer )
            {
                pCell = create_cell(*pCD);
                pCell->assign_position( x , y , 0.0 );

                if( fabs( y ) > 0.01 )
                {
                    pCell = create_cell(*pCD);
                    pCell->assign_position( x , -y , 0.0 );
                }

                if( fabs( x ) > 0.01 )
                {
                    pCell = create_cell(*pCD);
                    pCell->assign_position( -x , y , 0.0 );

                    if( fabs( y ) > 0.01 )
                    {
                        pCell = create_cell(*pCD);
                        pCell->assign_position( -x , -y , 0.0 );

                    }
                }
                x += cell_spacing;

            }

            y += cell_spacing * sqrt(3.0)/2.0;
            n++;
        }
	}

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

void CreateCell_InitialPhase( Cell* pCell, Phenotype& phenotype ){
	int idx_phase; // 0: G0G1 phase, 1: S phase, 2: G2 phase, 3: M phase, -1: Apoptotic
	#pragma omp critical
	{
		if ( InitialFraction[0] > 0 ){ // G0G1 phase
			idx_phase = 0;
			InitialFraction[0]--;
		}else if( InitialFraction[1] > 0 ){ // S phase
			idx_phase = 1;
			InitialFraction[1]--;
		}else if( InitialFraction[2] > 0 ){ // G2 phase
			idx_phase = 2;
			InitialFraction[2]--;
		}else if( InitialFraction[3] > 0 ){ // M phase
			idx_phase = 3;
			InitialFraction[3]--;
		}else{ // Apoptotic
			idx_phase = -1;
		}
		startingPhaseCytometricModel(pCell, phenotype, idx_phase );
	}
}

void cancer_phenotype( Cell* pCell, Phenotype& phenotype, double dt )
{
	if( phenotype.death.dead == true )
	{
		pCell->functions.update_phenotype = NULL;
		return;
	}
	if (dt == -1){ // Flag to generate initial condition with cells in different phases
		CreateCell_InitialPhase( pCell, phenotype);
		return;
	}
	int index_drug1 = microenvironment.find_density_index("drug1");
	int index_drug2 = microenvironment.find_density_index("drug2");
	pCell->custom_data["elapsed_time"] = pCell->phenotype.cycle.data.elapsed_time_in_phase;

	// Pressure negative feedback on proliferation
	double linear_response_pressure = decreasing_linear_response_function(get_single_signal(pCell, "pressure"), parameters.doubles("min_simple_pressure"), parameters.doubles("max_simple_pressure"));
	double entry_base_cycle_rate = get_single_base_behavior( pCell,"cycle entry");
	set_single_behavior( pCell,"cycle entry", entry_base_cycle_rate*linear_response_pressure);
	pCell->custom_data["simple_pressure"] = get_single_signal(pCell, "pressure");
	if ( parameters.bools["chemotherapy_drug1"].value ){
		// Pharmacokinetics of drug 1
		// Update total amount of drug inside of cell
		// std::cout << "old: " << pCell->custom_data["internal_drug1"] << " New: " << pCell->phenotype.molecular.internalized_total_substrates[index_drug1] << std::endl;
		pCell->custom_data["internal_drug1"] = pCell->phenotype.molecular.internalized_total_substrates[index_drug1];
		pCell->phenotype.molecular.internalized_total_substrates[index_drug1] -= dt*parameters.doubles("drug1_decay_rate")*pCell->custom_data["internal_drug1"];
		pCell->custom_data["internal_drug1"] = pCell->phenotype.molecular.internalized_total_substrates[index_drug1];

		// releasing efflux of drug 1
		double efflux_drug1 = parameters.doubles("drug1_efflux_rate")*pCell->custom_data["internal_drug1"]/pCell->phenotype.volume.total; // Export concentration not total amount (verified)
		set_single_behavior( pCell,"drug1 export", efflux_drug1);

		// Update drug
		double drug1_base_uptake = get_single_base_behavior(pCell, "drug1 uptake");
		double saturation_drug1_amount = pCell->phenotype.volume.total*parameters.doubles("internal_saturation_drug1"); // convert concentration to total amount
		double drug1_new_uptake =  drug1_base_uptake*decreasing_linear_response_function(pCell->custom_data["internal_drug1"], 0.1*saturation_drug1_amount,saturation_drug1_amount); // range of response is: 10% of base uptake to base uptake
		set_single_behavior( pCell,"drug1 uptake",drug1_new_uptake); // change the drug uptake according internalization

		// Saving new uptake and external drug1
		pCell->custom_data["uptake_drug1"] = drug1_new_uptake;
		double drug1_external_source = get_single_signal(pCell, "drug1"); // extracellular
		pCell->custom_data["external_drug1"] = drug1_external_source;

		// Pharmacodynamics of drug 1 (DNA damage-based mechanism)
		pCell->custom_data["dna_damage"] -= dt*parameters.doubles("dna_repair_rate")*pCell->custom_data["dna_damage"]; // repair dna damage in all cell cycle
		if ( pCell->phenotype.cycle.current_phase().code == PhysiCell_constants::S_phase ){ //check phase to damage damge (phase S)
			pCell->custom_data["dna_damage"] += dt*parameters.doubles("dna_damage_rate")*parameters.doubles("drug1_decay_rate")*pCell->custom_data["internal_drug1"]/pCell->phenotype.volume.total; // Based on concentration not in amount of internalized drug
		}

		// Apoptosis increase with DNA damage
		if ( pCell->custom_data["dna_damage"] > 0){
			double base_apoptosis = get_single_base_behavior( pCell, "apoptosis");
			double max_apoptosis = parameters.doubles("factor_max_apoptosis_drug1") * base_apoptosis;
			double hill_dna_damage = Hill_response_function( pCell->custom_data["dna_damage"] , parameters.doubles("dna_damage_apop_half_max") , parameters.doubles("dna_damage_hill_coef") ); // hill power = 1
			double apoptosis_rate = base_apoptosis + (max_apoptosis-base_apoptosis)*hill_dna_damage;
			set_single_behavior( pCell,"apoptosis",apoptosis_rate);
			pCell->custom_data["apoptosis_rate"] = apoptosis_rate;
		}
		// The Cell cycle arrest (end of phase M) increase with DNA damage and pressure
		if ( pCell->custom_data["dna_damage"] > 0 && pCell->phenotype.cycle.current_phase().code == PhysiCell_constants::M_phase && parameters.doubles("factor_arrest_MtoG0G1_drug1") > 0 ){
			double hill_dna_damage = Hill_response_function( pCell->custom_data["dna_damage"] , parameters.doubles("dna_damage_arrest_half_max") , parameters.doubles("dna_damage_hill_coef") );
			double base_M_G0G1_rate = get_single_base_behavior( pCell, "exit from cycle phase 3"); // transition phase M to G0/G1
			double fact_arrest_M_G0G1 = parameters.doubles("factor_arrest_MtoG0G1_drug1");
			double cycle_rate = base_M_G0G1_rate*(1-fact_arrest_M_G0G1*hill_dna_damage);
			set_single_behavior( pCell,"exit from cycle phase 3", cycle_rate);
			pCell->custom_data["rate_exit_M_phase"] = cycle_rate;
			// std::cout << "Drug 1 extra: " << get_single_signal(pCell, "drug1") << " drug 1 intra: " << drug1_internal_source << " total amont drug 1: " << pCell->custom_data["total_concentration_drug1"] << " Efflux: " << efflux_drug1 << std::endl;
			// std::cout << "DNA damage: " << pCell->custom_data["dna_damage"] << " cycle rate : " << cycle_rate << std::endl;
		}
	}

	if ( parameters.bools["chemotherapy_drug2"].value ){
		// Pharmacokinetics of drug 2
		// Update total amount of drug inside of cell
		pCell->phenotype.molecular.internalized_total_substrates[index_drug2] -= dt*parameters.doubles("drug2_decay_rate")*pCell->custom_data["internal_drug2"];
		pCell->custom_data["internal_drug2"] = pCell->phenotype.molecular.internalized_total_substrates[index_drug2];

		// releasing efflux of drug 2
		double efflux_drug2 = parameters.doubles("drug2_efflux_rate")*pCell->custom_data["internal_drug2"]/pCell->phenotype.volume.total; // Export concentration not total amount (verified);
		set_single_behavior( pCell,"drug2 export", efflux_drug2);

		// Update uptake
		double drug2_base_uptake = get_single_base_behavior(pCell, "drug2 uptake");
		double saturation_drug2_amount = pCell->phenotype.volume.total*parameters.doubles("internal_saturation_drug2"); // convert concentration to total amount
		double drug2_new_uptake =  drug2_base_uptake*decreasing_linear_response_function(pCell->custom_data["internal_drug2"], 0.1*saturation_drug2_amount,saturation_drug2_amount); // range of response is: 10% of base uptake to base uptake
		set_single_behavior( pCell,"drug2 uptake",drug2_new_uptake); // change the drug uptake according internalization

		// Saving new uptake and external drug1
		pCell->custom_data["uptake_drug2"] = drug2_new_uptake;
		double drug2_external_source = get_single_signal(pCell, "drug2"); // extracellular
		pCell->custom_data["external_drug2"] = drug2_external_source;

		// Pharmacodynamics of drug 2 (mitotic spindle-based mechanism)
		if ( pCell->phenotype.cycle.current_phase().code == PhysiCell_constants::M_phase ){  //check phase to damage and repair (phase M)
			pCell->custom_data["spindle_damage"] -= dt*parameters.doubles("spindle_repair_rate")*pCell->custom_data["spindle_damage"];
			pCell->custom_data["spindle_damage"] += dt*parameters.doubles("spindle_damage_rate")*parameters.doubles("drug2_decay_rate")*pCell->custom_data["internal_drug2"]/pCell->phenotype.volume.total; // Based on concentration not in amount of internalized drug;
			// Apoptosis increase with spindle damage
			double base_apoptosis = get_single_base_behavior( pCell, "apoptosis");
			double max_apoptosis = parameters.doubles("factor_max_apoptosis_drug2") * base_apoptosis;
			double hill_spindle_damage = Hill_response_function( pCell->custom_data["spindle_damage"] , parameters.doubles("spindle_damage_half_max") , parameters.doubles("spindle_damage_hill_coef") );
			double apoptosis_rate = base_apoptosis + (max_apoptosis-base_apoptosis)*hill_spindle_damage;
			set_single_behavior( pCell,"apoptosis",apoptosis_rate);
			pCell->custom_data["apoptosis_rate"] = apoptosis_rate;
		}//else don't change dna damage
		// The Cell cycle arrest (end of phase M) increase with spindle damage and pressure
		if ( pCell->custom_data["spindle_damage"] > 0 && pCell->phenotype.cycle.current_phase().code == PhysiCell_constants::M_phase){
			double hill_spindle_damage = Hill_response_function( pCell->custom_data["spindle_damage"] , parameters.doubles("spindle_damage_half_max") , parameters.doubles("spindle_damage_hill_coef") );
			double base_M_G0G1_rate = get_single_base_behavior( pCell, "exit from cycle phase 3"); // transition phase M to G0/G1
			double fact_arrest_M_G0G1 = parameters.doubles("factor_arrest_MtoG0G1_drug2");
			double cycle_rate = base_M_G0G1_rate*(1-fact_arrest_M_G0G1*hill_spindle_damage);
			set_single_behavior( pCell,"exit from cycle phase 3", cycle_rate);
			pCell->custom_data["rate_exit_M_phase"] = cycle_rate;
			// std::cout << "Drug 2 extra: " << get_single_signal(pCell, "drug2") << " drug 2 intra: " << drug2_internal_source << " total amont drug 2: " << pCell->custom_data["total_concentration_drug2"] << " Efflux: " << efflux_drug2 << std::endl;
			// std::cout << "Spindle damage: " << pCell->custom_data["spindle_damage"] << " cycle rate : " << cycle_rate << std::endl;
		}
	}

	return;
}

void create_plot_legend_custom( std::string filename, std::vector<std::string> (*cell_coloring_function)(Cell*) ){
	int number_of_cell_types = 5;

	double temp_cell_radius = 25;
	double temp_cell_volume = 4.1887902047863909846168578443727 * pow( temp_cell_radius , 3.0 );

	double relative_padding = 0.15;
	double padding = relative_padding * 2.0 * temp_cell_radius;

	double row_height = 2.0 * temp_cell_radius + 2*padding;

	double font_size = 0.85 * 2.0 * temp_cell_radius;
	double row_width  = 2.0 * temp_cell_radius + 2*padding + ( 32 * font_size ) + 2 * padding;

	double total_height = number_of_cell_types * row_height;
	double total_width  = row_width;

	std::ofstream os( filename , std::ios::out );
	Write_SVG_start( os , total_width ,total_height );

	double cursor_x = padding + temp_cell_radius;
	double cursor_y = padding + temp_cell_radius;

	static std::vector< std::string > output( 2 , "rgb(0,0,0)" );
	std::string name;
	for( int k=0 ; k < number_of_cell_types ; k++ )
	{
		if ( k == 0){
			name = "apoptotic";
			output[0] = "rgb(255,0,0)";
			output[1] = "rgb(125,0,0)";
		}
		if ( k == 1){
			name = "phase G0/G1";
			output[0] = "rgb(0,80,255)";
			output[1] = "rgb(0,40,255)";
		}
		if ( k == 2){
			name = "phase S";
			output[0] = "rgb(255, 0, 255)";
			output[1] = "rgb(190,0,190)";
		}
		if ( k == 3){
			name = "phase G2";
			output[0] = "rgb(255, 255, 0)";
			output[1] = "rgb(190, 190, 0)";
		}
		if ( k == 4){
		  name = "phase M";
			output[0] = "rgb(0,255,0)";
			output[1] = "rgb(0,190,0)";
		}

		// place a big circle with cytoplasm colors
		Write_SVG_circle(os,cursor_x, cursor_y , temp_cell_radius , 1.0 , "rgb(0,0,0)" , output[0] );
		// place a small circle with nuclear colors
		Write_SVG_circle(os,cursor_x, cursor_y , 0.5*temp_cell_radius , 1.0 , output[1] , "rgb(0,0,0)" );

		// place the label

		cursor_x += temp_cell_radius + 2*padding;
		cursor_y += 0.3*font_size;

		Write_SVG_text( os , name.c_str() , cursor_x , cursor_y, font_size ,
			PhysiCell_SVG_options.font_color.c_str() , PhysiCell_SVG_options.font.c_str() );

		// move the cursor down to the next row

		cursor_y -= 0.3*font_size;
		cursor_y += ( 2.0 * padding + 2.0*temp_cell_radius );
		cursor_x = padding + temp_cell_radius;
	}
	Write_SVG_end( os );
	os.close();
}

void CellCount(	std::vector<int>& NumberofCells, std::vector<double>& VolumeofCells, std::vector<double>& DNA_damage_quantiles )
{
	std::vector<double> DNA_damage_liveCells;
	for (int i=0; i < (*all_cells).size(); i++)
	{
		if( (*all_cells)[i]->phenotype.death.dead == true )
		{
			NumberofCells[4]++;
			VolumeofCells[4] += (*all_cells)[i]->phenotype.volume.total;
		}
		else if( (*all_cells)[i]->phenotype.cycle.current_phase().code == PhysiCell_constants::G0G1_phase )
		{
			NumberofCells[0]++;
			VolumeofCells[0] += (*all_cells)[i]->phenotype.volume.total;
			DNA_damage_liveCells.push_back((*all_cells)[i]->custom_data["dna_damage"]);
		}
		else if ( (*all_cells)[i]->phenotype.cycle.current_phase().code == PhysiCell_constants::S_phase )
		{
			NumberofCells[1]++;
			VolumeofCells[1] += (*all_cells)[i]->phenotype.volume.total;
			DNA_damage_liveCells.push_back((*all_cells)[i]->custom_data["dna_damage"]);
		}
		else if ( (*all_cells)[i]->phenotype.cycle.current_phase().code == PhysiCell_constants::G2_phase )
		{
			NumberofCells[2]++;
			VolumeofCells[2] += (*all_cells)[i]->phenotype.volume.total;
			DNA_damage_liveCells.push_back((*all_cells)[i]->custom_data["dna_damage"]);
		}
		else if ( (*all_cells)[i]->phenotype.cycle.current_phase().code == PhysiCell_constants::M_phase )
		{
			NumberofCells[3]++;
			VolumeofCells[3] += (*all_cells)[i]->phenotype.volume.total;
			DNA_damage_liveCells.push_back((*all_cells)[i]->custom_data["dna_damage"]);
		}
	}
	if (parameters.bools["chemotherapy_drug1"].value && DNA_damage_liveCells.size() > 10){
		// Sort vector in increasing order
		std::sort(DNA_damage_liveCells.begin(), DNA_damage_liveCells.end());
		int size = DNA_damage_liveCells.size();
		double q1_frac_part, q2_frac_part, q3_frac_part;
		double q1_int_part, q2_int_part, q3_int_part;
		q1_frac_part = modf (0.25*(size + 1) , &q1_int_part);
		q2_frac_part = modf (0.5*(size + 1) , &q2_int_part);
  	q3_frac_part = modf ((3.0/4.0)*(size + 1) , &q3_int_part);
  	float median, q1, q3;
		q1 = DNA_damage_liveCells[int(q1_int_part)] + q1_frac_part*(DNA_damage_liveCells[int(q1_int_part)+1]-DNA_damage_liveCells[int(q1_int_part)]);
  	median = DNA_damage_liveCells[int(q2_int_part)] + q2_frac_part*(DNA_damage_liveCells[int(q2_int_part)+1]-DNA_damage_liveCells[int(q2_int_part)]);
		q3 = DNA_damage_liveCells[int(q3_int_part)] + q3_frac_part*(DNA_damage_liveCells[int(q3_int_part)+1]-DNA_damage_liveCells[int(q3_int_part)]);
		DNA_damage_quantiles = {q1, median, q3};
		DNA_damage_liveCells.clear();
	}
}

void startingPhaseCytometricModel( Cell* pCell, Phenotype& phenotype, int idx_phase )
{
	if ( idx_phase > 0 ){
		pCell->phenotype.cycle.data.transition_rate(0,1)=9e99;
		pCell->phenotype.cycle.advance_cycle(pCell,phenotype, 0.01);
		set_single_behavior( pCell,"exit from cycle phase 0",get_single_base_behavior( pCell, "exit from cycle phase 0"));
		if ( idx_phase > 1 ){
			pCell->phenotype.cycle.data.transition_rate(1,2)=9e99;
			pCell->phenotype.cycle.advance_cycle(pCell,phenotype, 0.01);
			set_single_behavior( pCell,"exit from cycle phase 1",get_single_base_behavior( pCell, "exit from cycle phase 1"));
			if ( idx_phase == 3 ){
				pCell->phenotype.cycle.data.transition_rate(2,3)=9e99;
				pCell->phenotype.cycle.advance_cycle(pCell,phenotype, 0.01);
				set_single_behavior( pCell,"exit from cycle phase 2",get_single_base_behavior( pCell, "exit from cycle phase 2"));
			}
		}
	}
	if ( idx_phase == -1 ){
		static int apoptosis_index = phenotype.death.find_death_model_index( "Apoptosis" );
		pCell->start_death( apoptosis_index );
	}
	// std::cout << pCell->ID << ": " << pCell->phenotype.cycle.current_phase().name << " Dead: " << pCell->phenotype.death.dead << std::endl;
}

void InitializeCells( void ){
	// Set Fraction of each cell type
	InitialFraction[0] = (int)round( (*all_cells).size()*parameters.doubles("G0G1_initial_fraction") );
	InitialFraction[1] = (int)round( (*all_cells).size()*parameters.doubles("S_initial_fraction") );
	InitialFraction[2] = (int)round( (*all_cells).size()*parameters.doubles("G2_initial_fraction") );
	InitialFraction[3] = (int)round( (*all_cells).size()*parameters.doubles("M_initial_fraction") );
	// std::cout << InitialFraction[0] << " " << InitialFraction[1] << " " << InitialFraction[2] << " " << InitialFraction[3] << std::endl;
	// Update phase of each cell
	for( int i=0; i < (*all_cells).size(); i++ ){
		if( (*all_cells)[i]->is_out_of_domain == false ){
			(*all_cells)[i]->functions.update_phenotype( (*all_cells)[i] , (*all_cells)[i]->phenotype , -1.0 );
		}
	}
}

void Chemotherapy( std::string substrateName, double time_chemotherapy, double duration_chemotherapy, double dose_chemotherapy, double diffusion_dt ){
	// Check the time where the infusion is activated
	bool drug_activation;
	if( fabs( PhysiCell_globals.current_time - time_chemotherapy ) < 0.01 * diffusion_dt ){ // Add drug
		drug_activation = true;
		std::cout << "Starting chemotherapy with " << substrateName << " - Time: " << PhysiCell_globals.current_time << " Duration: " << duration_chemotherapy << " Dose: " << dose_chemotherapy <<std::endl;
	}
	else{
		if ( fabs( PhysiCell_globals.current_time - time_chemotherapy - duration_chemotherapy ) < 0.01 * diffusion_dt ){ // Stop drug releasing
			drug_activation = false;
			std::cout << "Finishing chemotherapy with " << substrateName << " - Time: " << PhysiCell_globals.current_time << std::endl;
		}
		else
			return; // It's not the period of the drug infusion
	}
	// Set the dose amount
	int subs_ind = microenvironment.find_density_index( substrateName );
	default_microenvironment_options.Dirichlet_xmin[subs_ind] = true; // Just in substrate 0 = drug 1
	default_microenvironment_options.Dirichlet_xmin_values[subs_ind] =  dose_chemotherapy;
	// Apply or Remove dirichlet node in x_min of domain - value = 1.0
	for( unsigned int k=0 ; k < microenvironment.mesh.z_coordinates.size() ; k++ ){
		int I = 0;
		// set Dirichlet conditions along the xmin outer edges
		for( unsigned int j=0 ; j < microenvironment.mesh.y_coordinates.size() ; j++ ){
			// set the value
			microenvironment.add_dirichlet_node( microenvironment.voxel_index(I,j,k) , default_microenvironment_options.Dirichlet_xmin_values );
			// set the activation or inactivation
			if (drug_activation){
				microenvironment.set_substrate_dirichlet_activation(  microenvironment.voxel_index(I,j,k) , default_microenvironment_options.Dirichlet_xmin );
			}
			else
				microenvironment.remove_dirichlet_node( microenvironment.voxel_index(I,j,k) );
		}
	}
}
