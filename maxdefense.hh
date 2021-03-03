///////////////////////////////////////////////////////////////////////////////
// maxdefense.hh
//
// Compute the set of armos that maximizes defense, within a gold budget,
// with the dynamic method or exhaustive search.
//
///////////////////////////////////////////////////////////////////////////////


#pragma once


#include <cassert>
#include <cmath>
#include <fstream>
#include <iostream>
#include <iomanip>
#include <memory>
#include <queue>
#include <sstream>
#include <string>
#include <vector>


// One armor item available for purchase.
class ArmorItem
{
	//
	public:
		
		//
		ArmorItem
		(
			const std::string& description,
			size_t cost_gold,
			double defense_points
		)
			:
			_description(description),
			_cost_gold(cost_gold),
			_defense_points(defense_points)
		{
			assert(!description.empty());
			assert(cost_gold > 0);
		}
		
		//
		const std::string& description() const { return _description; }
		int cost() const { return _cost_gold; }
		double defense() const { return _defense_points; }
	
	//
	private:
		
		// Human-readable description of the armor, e.g. "new enchanted helmet". Must be non-empty.
		std::string _description;
		
		// Cost, in units of gold; Must be positive
		int _cost_gold;
		
		// Defense points; most be non-negative.
		double _defense_points;
};


// Alias for a vector of shared pointers to ArmorItem objects.
typedef std::vector<std::shared_ptr<ArmorItem>> ArmorVector;


// Load all the valid armor items from the CSV database
// Armor items that are missing fields, or have invalid values, are skipped.
// Returns nullptr on I/O error.
std::unique_ptr<ArmorVector> load_armor_database(const std::string& path)
{
	std::unique_ptr<ArmorVector> failure(nullptr);
	
	std::ifstream f(path);
	if (!f)
	{
		std::cout << "Failed to load armor database; Cannot open file: " << path << std::endl;
		return failure;
	}
	
	std::unique_ptr<ArmorVector> result(new ArmorVector);
	
	size_t line_number = 0;
	for (std::string line; std::getline(f, line); )
	{
		line_number++;
		
		// First line is a header row
		if ( line_number == 1 )
		{
			continue;
		}
		
		std::vector<std::string> fields;
		std::stringstream ss(line);
		
		for (std::string field; std::getline(ss, field, '^'); )
		{
			fields.push_back(field);
		}
		
		if (fields.size() != 3)
		{
			std::cout
				<< "Failed to load armor database: Invalid field count at line " << line_number << "; Want 3 but got " << fields.size() << std::endl
				<< "Line: " << line << std::endl
				;
			return failure;
		}
		
		std::string
			descr_field = fields[0],
			cost_gold_field = fields[1],
			defense_points_field = fields[2]
			;
		
		auto parse_dbl = [](const std::string& field, double& output)
		{
			std::stringstream ss(field);
			if ( ! ss )
			{
				return false;
			}
			
			ss >> output;
			
			return true;
		};
		
		std::string description(descr_field);
		double cost_gold, defense_points;
		if (
			parse_dbl(cost_gold_field, cost_gold)
			&& parse_dbl(defense_points_field, defense_points)
		)
		{
			result->push_back(
				std::shared_ptr<ArmorItem>(
					new ArmorItem(
						description,
						cost_gold,
						defense_points
					)
				)
			);
		}
	}

	f.close();
	
	return result;
}


// Convenience function to compute the total cost and defense in an ArmorVector.
// Provide the ArmorVector as the first argument
// The next two arguments will return the cost and defense back to the caller.
void sum_armor_vector
(
	const ArmorVector& armors,
	int& total_cost,
	double& total_defense
)
{
	total_cost = total_defense = 0;
	for (auto& armor : armors)
	{
		total_cost += armor->cost();
		total_defense += armor->defense();
	}
}


// Convenience function to print out each ArmorItem in an ArmorVector,
// followed by the total kilocalories and protein in it.
void print_armor_vector(const ArmorVector& armors)
{
	std::cout << "*** Armor Vector ***" << std::endl;
	
	if ( armors.size() == 0 )
	{
		std::cout << "[empty armor list]" << std::endl;
	}
	else
	{
		for (auto& armor : armors)
		{
			std::cout
				<< "Ye olde " << armor->description()
				<< " ==> "
				<< "Cost of " << armor->cost() << " gold"
				<< "; Defense points = " << armor->defense()
				<< std::endl
				;
		}
		
		int total_cost;
		double total_defense;
		sum_armor_vector(armors, total_cost, total_defense);
		std::cout
			<< "> Grand total cost: " << total_cost << " gold" << std::endl
			<< "> Grand total defense: " << total_defense
			<< std::endl
			;
	}
}


// Convenience function to print out a 2D cache, composed of an std::vector<std::vector<double>>
// For sanity, will refuse to print a cache that is too large.
// Hint: When running this program, you can redirect stdout to a file,
//	which may be easier to view and inspect than a terminal
void print_2d_cache(const std::vector<std::vector<double>>& cache)
{
	std::cout << "*** 2D Cache ***" << std::endl;
	
	if ( cache.size() == 0 )
	{
		std::cout << "[empty]" << std::endl;
	}
	else if ( cache.size() > 250 || cache[1].size() > 250 )
	{
		std::cout << "[too large]" << std::endl;
	}
	else
	{
		for ( const std::vector<double> row : cache)
		{
			for ( double value : row )
			{
				std::cout << std::setw(5) << value;
			}
			std::cout << std::endl;
		}
	}
}

// Filter the vector source, i.e. create and return a new ArmorVector
// containing the subset of the armor items in source that match given
// criteria.
// This is intended to:
//	1) filter out armor with zero or negative defense that are irrelevant to our optimization
//	2) limit the size of inputs to the exhaustive search algorithm since it will probably be slow.
//
// Each armor item that is included must have at minimum min_defense and at most max_defense.
//	(i.e., each included armor item's defense must be between min_defense and max_defense (inclusive).
//
// In addition, the the vector includes only the first total_size armor items that match these criteria.
std::unique_ptr<ArmorVector> filter_armor_vector
(
	const ArmorVector& source,
	double min_defense,
	double max_defense,
	int total_size
)
{
	//1) filter out armor with zero or negative defense that are irrelevant to our optimization
	std::unique_ptr<ArmorVector> result(new ArmorVector);
	int count = 0;
	for (auto& armor : source) {
		if (count < total_size) {
			if (armor->defense() > 0) {
				if (armor->defense() >= min_defense && armor->defense() <= max_defense) {
					result->push_back(armor);
					count++;
				}
			}
		}
		else {
			break;
		}
	}
	return result;
}

// To heapify a subtree rooted with node i which is 
// an index in arr[]. n is size of heap 
void heapify(ArmorVector& todo, int n, int i)
{
	int smallest = i; // Initialize largest as root 
	int l = 2 * i + 1; // left = 2*i + 1 
	int r = 2 * i + 2; // right = 2*i + 2 

	// If left child is larger than root 
	if (l < n && todo[l]->defense() / todo[l]->cost() < todo[smallest]->defense() / todo[smallest]->cost())
		smallest = l;

	// If right child is larger than largest so far 
	if (r < n && todo[r]->defense() / todo[r]->cost() < todo[smallest]->defense() / todo[smallest]->cost())
		smallest = r;

	// If largest is not root 
	if (smallest != i) {
		auto temp = todo[i];
		todo[i] = todo[smallest];
		todo[smallest] = temp;

		// Recursively heapify the affected sub-tree 
		heapify(todo, n, smallest);
	}
}

// main function to do heap sort 
void heapSort(ArmorVector& todo, int n)
{
	// Build heap (rearrange array) 
	for (int i = n / 2 - 1; i >= 0; i--) {
		heapify(todo, n, i);
	}

	// One by one extract an element from heap 
	for (int i = n - 1; i >= 0; i--) {
		// Move current root to end 
		auto temp = todo[0];
		todo[0] = todo[i];
		todo[i] = temp;

		// call max heapify on the reduced heap 
		heapify(todo, i, 0);
	}
}

// Define structure of a cell of a Dynamic Programming Table
struct dpstruct {
	double dpdefense;
	std::vector< int > arr;
};

// Compute the optimal set of armor items with a dynamic algorithm.
// Specifically, among the armor items that fit within a total_cost gold budget,
// choose the selection of armors whose defense is greatest.
// Repeat until no more armor items can be chosen, either because we've run out of armor items,
// or run out of gold.
std::unique_ptr<ArmorVector> dynamic_max_defense
(
	const ArmorVector& armors,
	int total_cost
)
{
	//todo = armor_items
	std::unique_ptr<ArmorVector> todo(new ArmorVector);
	for (auto& armor : armors) {
		todo->push_back(armor);
	}
	int m = total_cost;
	int n = todo->size();

	//heap sort
	heapSort(*todo, n);

	std::unique_ptr<ArmorVector> finalresult(new ArmorVector);

	//allocate the array
	dpstruct dp[3][m + 1];

	//set default
	dp[0][0].dpdefense = 0.0;
	for (int j = 1; j <= total_cost; j++) {
		dp[0][j].dpdefense = 0.0;
	}
	for (int i = 1; i <= 2; i++) {
		dp[i][0].dpdefense = 0.0;
	}

	int i = 1;
	int newi, previ;
	newi = i;
	previ = i - 1;

	//make dynamic knapsack table
	while (i <= n) {
		int itemCost = (*todo)[i - 1]->cost();
		double itemDefense = (*todo)[i - 1]->defense();

		for (int j = 1; j <= total_cost; j++) {
			//if itemCost > j, dp[i][j] = dp[i-1][j] 
			if (itemCost > j) {
				dp[newi][j].dpdefense = dp[previ][j].dpdefense;
				dp[newi][j].arr = dp[previ][j].arr;
			}
			else {
				//if dp[i][j] = max(dp[i-1][j].defense, itemDefense + dp[i-1][j - itemCost].dpdefense, add into dp[i][j]
				if (dp[previ][j].dpdefense > itemDefense + dp[previ][j - itemCost].dpdefense) {
					dp[newi][j].dpdefense = dp[previ][j].dpdefense;
					dp[newi][j].arr = dp[previ][j].arr;
				}
				else {
					dp[newi][j].dpdefense = itemDefense + dp[previ][j - itemCost].dpdefense;
					dp[newi][j].arr = dp[previ][j - itemCost].arr;
					dp[newi][j].arr.push_back(i - 1);
				}
			}
			//add result into finalresult
			if ((i == n) && (j == total_cost)) {
				std::vector< int > arr2;
				arr2 = dp[newi][j].arr;
				while (arr2.size() != 0) {
					int r = (*arr2.begin());
					finalresult->insert(finalresult->begin(), (*todo)[r]);
					arr2.erase(arr2.begin());
				}
			}
		}

		//switch new and previous space, clear privious space
		if ((previ == 0) || (previ == 2))
		{
			newi = 2;
			previ = 1;
		}
		else if (previ == 1)
		{
			newi = 1;
			previ = 2;
		}
		else
		{
			//std::cout << "\n463 newi =" << newi << ", previ = " << previ << "\n";
		} 

		for (int t = 1; t <= total_cost; t++) {
			dp[newi][t].dpdefense = 0.0;
			dp[newi][t].arr.clear();
		}

		i++;
	}

	//return finalresult
	return finalresult;
}


// Compute the optimal set of armor items with an exhaustive search algorithm.
// Specifically, among all subsets of armor items,
// return the subset whose gold cost fits within the total_cost budget,
// and whose total defense is greatest.
// To avoid overflow, the size of the armor items vector must be less than 64.
std::unique_ptr<ArmorVector> exhaustive_max_defense
(
	const ArmorVector& armors,
	double total_cost
)
{
	//n = |armor_items|
	const int n = armors.size();
	assert(n < 64);

	//best = None
	std::unique_ptr<ArmorVector> best(new ArmorVector);
	double best_defense = 0.0;

	//for bits from 0 to (2^n -1):
	for (int bits = 0;bits < pow(2, n);bits++) {
		//candidate = empty vector
		std::unique_ptr<ArmorVector> candidate(new ArmorVector);
		double candidate_total_cost = 0.0;
		double candidate_total_defense = 0.0;
		//for j from 0 to n-1:
		for (int j = 0;j < n;j++) {
			//if ((bits >> j) & 1) == 1:
			if (((bits >> j) & 1) == 1) {
				//candidate.add_back(armor_items[j])
				candidate->push_back(armors[j]);
				candidate_total_cost += armors[j]->cost();
				candidate_total_defense += armors[j]->defense();
			}
		}
		//if total_gold_cost(candidate) <= G:
		if (candidate_total_cost <= total_cost) {
			//if best is None or total_defense(candidate) > total_defense(best) :
			if (best->size() == 0 || candidate_total_defense > best_defense) {
				best->clear();
				//best = candidate;
				for (auto& armor : *candidate) {
					best->push_back(armor);
					best_defense = candidate_total_defense;
				}
			}
		}
	}
	return best;
}









