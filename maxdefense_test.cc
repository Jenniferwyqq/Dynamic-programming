///////////////////////////////////////////////////////////////////////////////
// maxdefense_test.cc
//
// Unit tests for maxdefense.hh
//
///////////////////////////////////////////////////////////////////////////////


#include <cassert>
#include <sstream>


#include "maxdefense.hh"
#include "rubrictest.hh"


int main()
{
	Rubric rubric;
	
	ArmorVector trivial_armors;
	trivial_armors.push_back(std::shared_ptr<ArmorItem>(new ArmorItem("test helmet", 10.0, 20.0)));
	trivial_armors.push_back(std::shared_ptr<ArmorItem>(new ArmorItem("test boots", 4.0, 5.0)));
	
	auto all_armors = load_armor_database("armor.csv");
	assert( all_armors );
	
	auto filtered_armors = filter_armor_vector(*all_armors, 1, 2500, all_armors->size());
	
	//
	rubric.criterion(
		"load_armor_database still works", 2,
		[&]()
		{
			TEST_TRUE("non-null", all_armors);
			TEST_EQUAL("size", 8064, all_armors->size());
		}
	);
	
	//
	rubric.criterion(
		"filter_armor_vector", 2,
		[&]()
		{
			auto
				three = filter_armor_vector(*all_armors, 100, 500, 3),
				ten = filter_armor_vector(*all_armors, 100, 500, 10);
			TEST_TRUE("non-null", three);
			TEST_TRUE("non-null", ten);
			TEST_EQUAL("total_size", 3, three->size());
			TEST_EQUAL("total_size", 10, ten->size());
			TEST_EQUAL("contents", "used high-quality mystical human chest plate", (*ten)[0]->description());
			TEST_EQUAL("contents", "deteriorating poor quality enchanted human shield", (*ten)[9]->description());
			for (int i = 0; i < 3; i++) {
				TEST_EQUAL("contents", (*three)[i]->description(), (*ten)[i]->description());
			}
		}
	);
	
	//
	rubric.criterion(
		"dynamic_max_defense trivial cases", 2,
		[&]()
		{
			std::unique_ptr<ArmorVector> soln;
			
			soln = dynamic_max_defense(trivial_armors, 3);
			TEST_TRUE("non-null", soln);
			TEST_TRUE("empty solution", soln->empty());
			
			soln = dynamic_max_defense(trivial_armors, 10);
			TEST_TRUE("non-null", soln);
			TEST_EQUAL("helmet only", 1, soln->size());
			TEST_EQUAL("helmet only", "test helmet", (*soln)[0]->description());
			
			soln = dynamic_max_defense(trivial_armors, 9);
			TEST_TRUE("non-null", soln);
			TEST_EQUAL("boots only", 1, soln->size());
			TEST_EQUAL("boots only", "test boots", (*soln)[0]->description());
			
			soln = dynamic_max_defense(trivial_armors, 14);
			TEST_TRUE("non-null", soln);
			TEST_EQUAL("helmet and boots", 2, soln->size());
			TEST_EQUAL("helmet and boots", "test boots", (*soln)[0]->description());
			TEST_EQUAL("helmet and boots", "test helmet", (*soln)[1]->description());
		}
	);
	
	//
	rubric.criterion(
		"dynamic_max_defense correctness", 4,
		[&]()
		{
			std::unique_ptr<ArmorVector> soln_small, soln_large;
			
			soln_small = dynamic_max_defense(*filtered_armors, 500),
			soln_large = dynamic_max_defense(*filtered_armors, 5000);
			
			//print_armor_vector(*soln_small);
			//print_armor_vector(*soln_large);
			
			TEST_TRUE("non-null", soln_small);
			TEST_TRUE("non-null", soln_large);
			
			TEST_FALSE("non-empty", soln_small->empty());
			TEST_FALSE("non-empty", soln_large->empty());
			
			int
				cost_small,
				cost_large
				;
			double
				defense_small,
				defense_large
				;
			sum_armor_vector(*soln_small, cost_small, defense_small);
			sum_armor_vector(*soln_large, cost_large, defense_large);
			
			//	Precision
			defense_small	= std::round( defense_small	* 100 ) / 100;
			defense_large	= std::round( defense_large	* 100 ) / 100;
			
			TEST_EQUAL("Small solution cost",	500,		cost_small);
			TEST_EQUAL("Small solution defense",	9564.92,	defense_small);
			TEST_EQUAL("Large solution cost",	5000,		cost_large);
			TEST_EQUAL("Large solution defense",	82766.449999999997,	defense_large);
		}
	);
	
	//
	rubric.criterion(
		"exhaustive_max_defense trivial cases", 2,
		[&]()
		{
			std::unique_ptr<ArmorVector> soln;
			
			soln = exhaustive_max_defense(trivial_armors, 3);
			TEST_TRUE("non-null", soln);
			TEST_TRUE("empty solution", soln->empty());
			
			soln = exhaustive_max_defense(trivial_armors, 10);
			TEST_TRUE("non-null", soln);
			TEST_EQUAL("helmet only", 1, soln->size());
			TEST_EQUAL("helmet only", "test helmet", (*soln)[0]->description());
			
			soln = exhaustive_max_defense(trivial_armors, 9);
			TEST_TRUE("non-null", soln);
			TEST_EQUAL("boots only", 1, soln->size());
			TEST_EQUAL("boots only", "test boots", (*soln)[0]->description());
			
			soln = exhaustive_max_defense(trivial_armors, 14);
			TEST_TRUE("non-null", soln);
			TEST_EQUAL("helmet and boots", 2, soln->size());
			TEST_EQUAL("helmet and boots", "test helmet", (*soln)[0]->description());
			TEST_EQUAL("helmet and boots", "test boots", (*soln)[1]->description());
		}
	);
	
	//
	rubric.criterion(
		"exhaustive_max_defense correctness", 4,
		[&]()
		{
			std::vector<double> optimal_defense_totals =
			{
				500,		1033.05,	1500,	2100,	2400,
				2900,		3400,		4200,	4300,	4600,
				5000,		5400,		5800,	6100,	6500,
				7000,		7500,		8100,	8600,	8700
			};
			
			for ( size_t optimal_index = 0; optimal_index < optimal_defense_totals.size(); optimal_index++ )
			{
				size_t n = optimal_index + 1;
				double expected_defense = optimal_defense_totals[optimal_index];
				
				auto small_armors = filter_armor_vector(*filtered_armors, 1, 2000, n);
				TEST_TRUE("non-null", small_armors);
				
				auto solution = exhaustive_max_defense(*small_armors, 2000);
				TEST_TRUE("non-null", solution);
				
				int actual_cost;
				double actual_defense;
				sum_armor_vector(*solution, actual_cost, actual_defense);
				
				// Round
				expected_defense	= std::round( expected_defense	/ 100.0) * 100;
				actual_defense		= std::round( actual_defense	/ 100.0) * 100;
				
				std::stringstream ss;
				ss
					<< "exhaustive search n = " << n << " (optimal index = " << optimal_index << ")"
					<< ", expected defense = " << expected_defense
					<< " but algorithm found = " << actual_defense
					;
				TEST_EQUAL(ss.str(), expected_defense, actual_defense);
				
				auto dynamic_solution = dynamic_max_defense(*small_armors, 2000);
				int dynamic_actual_cost;
				double dynamic_actual_defense;
				sum_armor_vector(*solution, dynamic_actual_cost, dynamic_actual_defense);
				dynamic_actual_defense	= std::round( dynamic_actual_defense	/ 100.0) * 100;
				TEST_EQUAL("Exhaustive and dynamic get the same answer", actual_defense, dynamic_actual_defense);
			}
		}
	);
	
	return rubric.run();
}




