// compile & run:     g++ main.cpp -o main && ./main > out.txt
// © Jose Garfias Lopez
// Genetic  Algorithm

#include <iostream>
#include <cstdio>
#include <vector>
#include <math.h>

using namespace std; 

#define EXPERIMENT_NUMBER 3
#define GENERATIONS_TO_RUN 3000
#define POPULATION_SIZE 100
#define SIZE_OF_CHROMOSOME 390
#define SIZE_OF_CHUNK 13
const string GENES = "01"; 


///////////////////////////   HELPER FUNCTIONS   ///////////////////////////////////////////

int random_num(int start, int end) {
    int range = (end - start) + 1; 
    int random_int = start + (rand() % range); 
    return random_int; 
}

float random_float_num(float start, float end) {
    int range = (end - start) + 1; 
    float random_float = start + (rand() % range); 
    return random_float; 
}

int binary_to_decimal(string n) {
    string num = n; 
    int dec_value = 0; 
    int base = 1;  
    int len = num.length(); 
    for (int i = len - 1; i >= 0; i--) {
        if (num[i] == '1') {
            dec_value += base; 
        }
        base = base * 2; 
    }
    return dec_value; 
}

char mutated_genes() { 
    int r = random_num(0, GENES.size() - 1); 
    return GENES[r];
} 

string create_gnome() {
    string gnome = "";
    for(int i=0; i < SIZE_OF_CHROMOSOME; i++) {
        gnome += mutated_genes();
    }
    return gnome; 
}

float chromosome_to_x(string chromosome) {
    float aj = -30;
    float bj = 30;
    return aj + (float)binary_to_decimal(chromosome) * ((bj - aj) / ((pow(2,13) - 1)));
}
//////////////////////////////////////////////////////////////////////
  
// Class representing individual in population 
class Individual {
public: 
    string chromosome; 
    vector<float> x_values;
    float fitness;
    float p_select;
    float expected_value;
    int id;
    int generation;

    Individual(string chromosome, int id, int generation); 
    vector<Individual> mate(Individual parent2, int id1, int id2); 
    void calc_fitness(); 
    void calc_x_values(); 
    void uniform_mutate(); 
};
  
Individual::Individual(string chromosome, int id, int generation) {
    this->chromosome = chromosome; 
    this->calc_fitness();  
    this->id = id;
    this->generation = generation;
}; 

vector<Individual> Individual::mate(Individual partner, int id1, int id2) {
    vector<Individual> childs;
    string parentA_A = "";
    string parentA_B = "";
    string parentB_A = "";
    string parentB_B = "";
    int point_of_cross = random_num(1, SIZE_OF_CHROMOSOME - 1);

    for (int i=0; i <= point_of_cross; i++) {
        parentA_A += this->chromosome[i];
        parentB_A += partner.chromosome[i];
    }
    for (int i=point_of_cross+1; i < SIZE_OF_CHROMOSOME; i++) {
        parentA_B += this->chromosome[i];
        parentB_B += partner.chromosome[i];
    }

    childs.push_back(Individual(parentA_A + parentB_B, id1, this->generation + 1));
    childs.push_back(Individual(parentA_B + parentB_A, id2, this->generation + 1));

    return childs; 
};

void Individual::calc_x_values () {
    // CALCULATE X PER SUBCHROMOSOME IN CHROMOSOME
    this->x_values.clear();
    int chunk = 0;
    string sub_chromosome = "";
    for (int i=0; i<SIZE_OF_CHROMOSOME; i++) {
        sub_chromosome += this->chromosome[i];
        if (sub_chromosome.length() == 13) {
            this->x_values.push_back(chromosome_to_x(sub_chromosome));
            sub_chromosome = "";
        }
    }
}

void Individual::calc_fitness() {
    this->calc_x_values();

    float sum_result = 0;
    for (int i=0; i<this->x_values.size(); i++) {
        float Xi = x_values[i];
        sum_result += abs(100 * powf((Xi+1)-(powf(Xi, 2)), 2) + powf((Xi-1), 2));
    }

    this->fitness = sum_result;
};

void Individual::uniform_mutate() {
    int random_position = random_num(0, SIZE_OF_CHROMOSOME - 1);
    if (this->chromosome[random_position] == '0') {
        this->chromosome[random_position] = '1';
    } else {
        this->chromosome[random_position] = '0';
    }
    this->calc_fitness();
};

bool operator<(const Individual &ind1, const Individual &ind2) {
    return ind1.fitness > ind2.fitness; 
}

void calc_expected_values(vector<Individual> &population) {
    int n = population.size();
    float sumOfFitness = 0;
    for (int i=0; i<n; i++) {
        sumOfFitness += population[i].fitness;
    }
    for (int i=0; i<n; i++) {
        population[i].p_select = ((population[i].fitness) / sumOfFitness);
        population[i].expected_value = population[i].p_select * (float)n;
    }
}

int roulete_selection (vector<Individual> &population) {
    int n = population.size();
    float T = 0;
    for (int i=0; i<n; i++) {
        T += population[i].expected_value;
    }

    float r = random_float_num(0.0, T);
    int individual_i = -1;
    float sum=0;
    
    for (int i=0; i<n; i++) {
        sum += population[i].expected_value;
        if (sum >= r) {
            individual_i = i;
            break;
        }
    }

    return individual_i; 
}

void print_headers() {
    if (EXPERIMENT_NUMBER == 1) {
        cout<<"generacion, id, cromosoma, ";
        for (int i=0; i<(SIZE_OF_CHROMOSOME/SIZE_OF_CHUNK); i++) {
            cout << "valor " << i+1 <<",";
        }
        cout<<"aptitud, valor_esperado,"<<endl;
    }

    if (EXPERIMENT_NUMBER == 2) {
        cout<<"generacion, id, cromosoma, ";
        for (int i=0; i<(SIZE_OF_CHROMOSOME/SIZE_OF_CHUNK); i++) {
            cout << "valor " << i+1 <<",";
        }
        cout<<"aptitud_individuo, valor_esperado, promedio_generacion"<<endl;
    }

    if (EXPERIMENT_NUMBER == 3) {
        cout<<"generacion, id, cromosoma, ";
        for (int i=0; i<(SIZE_OF_CHROMOSOME/SIZE_OF_CHUNK); i++) {
            cout << "valor " << i+1 <<",";
        }
        cout<<"aptitud_individuo, valor_esperado, promedio_generacion"<<endl;
    }
}

void print_generation(int genration_number, vector<Individual> &population) {
    if (EXPERIMENT_NUMBER == 1 && genration_number != 0) {
        for (int i=0; i<population.size(); i++) {
            cout<< genration_number << ","; 
            cout<< population[i].id << ","; 
            cout<< population[i].chromosome <<","; 
            for (int j=0; j<(SIZE_OF_CHROMOSOME/SIZE_OF_CHUNK); j++) {
                cout << population[i].x_values[j] <<",";
            }
            cout<< population[i].fitness << ",";
            cout<< population[i].expected_value << ",\n";
        }
        cout<<",,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,\n";
        cout<<",,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,\n";
        cout<<",Mejor candidato,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,\n";
        cout<< genration_number << ","; 
        cout<< population[0].id << ","; 
        cout<< population[0].chromosome <<","; 
        for (int j=0; j<(SIZE_OF_CHROMOSOME/SIZE_OF_CHUNK); j++) {
            cout << population[0].x_values[j] <<",";
        }
        cout<< population[0].fitness << ",";
        cout<< population[0].expected_value << ",\n";
        cout<<",,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,\n";
        cout<<",,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,\n";
        cout<<",Padres seleccionados,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,\n";
    }

    if (EXPERIMENT_NUMBER == 2) {
        cout<< genration_number << ","; 
        cout<< population[0].id << ","; 
        cout<< population[0].chromosome <<","; 
        for (int j=0; j<(SIZE_OF_CHROMOSOME/SIZE_OF_CHUNK); j++) {
            cout << population[0].x_values[j] <<",";
        }
        cout<< population[0].fitness << ",";
        cout<< population[0].expected_value << ",";

        float sum = 0;
        for (int i=0; i<population.size(); i++) {
            sum+=population[i].fitness;
        }
        cout<< sum/population.size() << ",\n";
    }

    if (EXPERIMENT_NUMBER == 3) {
        cout<< genration_number << ","; 
        cout<< population[0].id << ","; 
        cout<< population[0].chromosome <<","; 
        for (int j=0; j<(SIZE_OF_CHROMOSOME/SIZE_OF_CHUNK); j++) {
            cout << population[0].x_values[j] <<",";
        }
        cout<< population[0].fitness << ",";
        cout<< population[0].expected_value << ",";

        float sum = 0;
        for (int i=0; i<population.size(); i++) {
            sum+=population[i].fitness;
        }
        cout<< sum/population.size() << ",\n";
    }
}
void print_parents(int genration_number, Individual p1, Individual p2) {
    if (EXPERIMENT_NUMBER == 1 && genration_number != 0) {
        cout<< genration_number << ","; 
        cout<< p1.id << ","; 
        cout<< p1.chromosome <<","; 
        for (int j=0; j<(SIZE_OF_CHROMOSOME/SIZE_OF_CHUNK); j++) {
            cout << p1.x_values[j] <<",";
        }
        cout<< p1.fitness << ",";
        cout<< p1.expected_value << ",\n";
        
        cout<< genration_number << ","; 
        cout<< p2.id << ","; 
        cout<< p2.chromosome <<","; 
        for (int j=0; j<(SIZE_OF_CHROMOSOME/SIZE_OF_CHUNK); j++) {
            cout << p2.x_values[j] <<",";
        }
        cout<< p2.fitness << ",";
        cout<< p2.expected_value << ",\n";

        cout<<",,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,\n";
    }
}
void print_better_worst_of_all_gens(Individual better, Individual worst) {
    if (EXPERIMENT_NUMBER == 3) {
        cout<<",,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,\n";
        cout<<",,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,\n";

        cout<<",,Mejor de todas las generaciones,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,\n";
        cout<< better.generation << ","; 
        cout<< better.id << ","; 
        cout<< better.chromosome <<","; 
        for (int j=0; j<(SIZE_OF_CHROMOSOME/SIZE_OF_CHUNK); j++) {
            cout << better.x_values[j] <<",";
        }
        cout<< better.fitness << ",";
        cout<< ",,\n";
        cout<<",,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,\n";
        cout<<",,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,\n";


        cout<<",,PEOR de todas las generaciones,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,\n";
        cout<< worst.generation << ","; 
        cout<< worst.id << ","; 
        cout<< worst.chromosome <<","; 
        for (int j=0; j<(SIZE_OF_CHROMOSOME/SIZE_OF_CHUNK); j++) {
            cout << worst.x_values[j] <<",";
        }
        cout<< worst.fitness << ",";
        cout<<",,\n";
        cout<<",,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,\n";
        cout<<",,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,\n";
    }
}
  
int main() {
    print_headers();
    srand((unsigned)(time(0)));
    int generation = 0; 
    vector<Individual> population; 
    vector<Individual> new_generation;

    // create initial population
    for(int i=0; i<POPULATION_SIZE; i++) {
        string gnome = create_gnome();
        population.push_back(Individual(gnome, i, generation));
    }
    calc_expected_values(population);
    sort(population.begin(), population.end());

    Individual betterIndividual = Individual(population[0].chromosome, population[0].id, population[0].generation);
    Individual worstIndividual = Individual(population[0].chromosome, population[0].id, population[0].generation);

    while(generation <= GENERATIONS_TO_RUN) {
        // SORT THE GENERATION BY FITNESS
        sort(population.begin(), population.end());

        if (population[0].fitness > betterIndividual.fitness) {
            betterIndividual = Individual(population[0].chromosome, population[0].id, population[0].generation);
        }
        if (population[0].fitness < worstIndividual.fitness) {
            worstIndividual = Individual(population[0].chromosome, population[0].id, population[0].generation);
        }

        print_generation(generation, population);

        new_generation.clear(); // clear the array for new generation

        // ELITISM ON THE FIRST 20% BEST CANDIDATES
        int s = (20*POPULATION_SIZE)/100; 
        for(int i = 0; i<s; i++){
            new_generation.push_back(population[i]);
        }
        vector<Individual> remaining_population;
        for(int i=s; i<population.size(); i++){
            remaining_population.push_back(population[i]);
        }
    
        // SELECTION AND CROSSOVER OVER THE REMAINING 80% 
        s = ((80*POPULATION_SIZE)/100 ) / 2;
        for(int i = 0; i<s; i++) {
            Individual parent1 = remaining_population[roulete_selection(remaining_population)];
            Individual parent2 = remaining_population[roulete_selection(remaining_population)];
            print_parents(generation, parent1, parent2);
            vector<Individual> childs = parent1.mate(parent2, new_generation.size(), new_generation.size() +1);
            new_generation.push_back(childs[0]);
            new_generation.push_back(childs[1]);
        }

        // MUTATION ON ALL: 1/L
        for (int i=0; i<new_generation.size(); i++) {
            new_generation[i].uniform_mutate();
        }

        population = new_generation;
        calc_expected_values(population);
        generation++;
    }
    
    print_better_worst_of_all_gens(betterIndividual, worstIndividual);
}