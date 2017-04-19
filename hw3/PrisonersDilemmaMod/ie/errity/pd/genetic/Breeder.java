package ie.errity.pd.genetic;

import ie.errity.pd.*;
import javax.swing.*;
import java.awt.*;
import java.awt.event.*;
import java.util.BitSet;
import java.util.Random;
import java.util.Collections;
import java.util.Arrays;


/**
 * Provides a means of evolving a population of 
 * {@link  ie.errity.pd.Prisoner Prisoners}
 * via a genetic algorithm
 * @author Andrew Errity 99086921
 * @author Sherri Goings (modified 4/4/2017)
 * @author Hazel Que
 */
public class Breeder extends JPanel
{
    private Prisoner curPopulation[];
    private double mutateP, crossP; //mutate, cross probabilities
    private int selection, selParam, seed; // selection method, any associated selection parameter, random number seed
    private int popSize;
    private Random rand;
	
    /**
     *Create a new Genetic Breeder
     */
    public Breeder()
    {	
		//defaults
		mutateP = .001;
		crossP = .95;
		selection = 0;
		selParam = 1;
		seed = -1;
		rand = new Random(); //uses current time as seed
    }
	
    /**
     *Create a new Genetic Breeder using {@link  ie.errity.pd.Rules Rules} given
     *@param r1 parameters for the breeder
     */
    public Breeder(Rules r1)
    {	
		//init based on rules
		mutateP = r1.getMutateP();
		crossP = r1.getCrossP();
		selection = r1.getSelection();
		selParam = r1.getSelectionParameter();
		seed = r1.getSeed();
		if (seed==-1)
			rand = new Random();
		else
			rand = new Random(seed); //uses seed set in rules for repeatability
    }
	
    /**
     *Breeds the next generation (panmictic mating) of an array of 
     *{@link  ie.errity.pd.Prisoner Prisoners}
     *@param c	initial population (raw fitness of population must be calculated previously)
     *@return the next generation
     */
    public Prisoner[] Breed(Prisoner[] c) {
		curPopulation = c;	//population to breed
		popSize = curPopulation.length;
		Prisoner Selected[] = new Prisoner[popSize]; // parent pop after selection


		// Select parents for next gen
		// ***ADD CODE (make separate functions) to perform each of the types of selection***
		// selection is an int that determines the method to use
		// selParam is an int that can be used as a parameter for a selection method if required

		// One silly selection method which uses rand and the selParam is given below
		// Be sure to use "rand" class variable here as your random number generator (i.e. don't create a new one!)

		// selection method 0, use parameter as a threshhold percentage.  If a random number is less than param,
		// select the best individual from the population.  Otherwise select a random individual.
		// In this method a param of 0 gives random selection, and a param of 100 gives best-wins-all selection.
		if (selection == 0) {
			// find index of most fit individual
			double maxFit = 0;
			int indexBest = 0;
			for (int i=0; i<popSize; i++) {
				if (curPopulation[i].getScore() > maxFit) {
					maxFit = curPopulation[i].getScore();
					indexBest = i;
				}
			}

			// select according to description above for this method
			for (int i=0; i<popSize; i++) {
			int selIndex = 0;
			if (rand.nextInt(100) < selParam)  // rand less than param, select best individual
				{
				selIndex = indexBest;
				}
			else  // otherwise select random individual
				{
				selIndex = rand.nextInt(popSize);
				}
			Selected[i] = (Prisoner)curPopulation[selIndex].clone();
			}
		}

		/* fitness proportional selection with sigma scaling and stochastic
		universal sampling, plus optional elitism. */
		else if(selection == 1) {
			Selected = fitnessProportionateSelection(curPopulation, selParam);
		}

		/* tournament selection */
		else if(selection == 2) {
			Selected = tournamentSelection(curPopulation, selParam);
		}

		else {  // any other selection method fill pop with always cooperate
			for (int i=0; i<popSize; i++)
			Selected[i] = new Prisoner("ALLC");
		}


		//Crossover & Mutate each pair of selected parents
		BitSet Offspring[] = new BitSet[2];  // temporarily holds 2 children during crossover/mutation
		for (int d=0; d<popSize; d+=2) {
			// in case of odd population, just mutate and replace last individual
			if (d+1 >= popSize) {
				Offspring[0] = Genetic.mutate(Selected[d].getStrat(), mutateP, rand);
				Selected[d] = new Prisoner(Offspring[0]);
			}
			else {

				if(rand.nextDouble() <= crossP) { //Cross Over
					Offspring = Genetic.crossover(Selected[d].getStrat(), Selected[d + 1].getStrat(), rand);
				}

				else  { //clones
					Offspring[0] = (BitSet)Selected[d].getStrat().clone();
					Offspring[1] = (BitSet)Selected[d+1].getStrat().clone();
				}

				//Mutation
				Offspring[0] = Genetic.mutate(Offspring[0],mutateP, rand);
				Offspring[1] = Genetic.mutate(Offspring[1],mutateP, rand);

				// Replacement - we are done with parents d & d+1, so just replace with children without
				// creating an entire new array
				Selected[d] = new Prisoner(Offspring[0]);
				Selected[d+1] = new Prisoner(Offspring[1]);
			}
		}
		// pass on children pop to be parents of next gen
		curPopulation = Selected;
		repaint();	//update display (if any)
		return curPopulation; //return the bred population
    }

	/**
	 * Enables fitness proportional selection with sigma scaling and stochastic
	 * universal sampling, plus optional elitism
	 * @param c - initial population
	 * @param m - elitism disabled if m is 0. otherwise, number of most fit
	 *          individuals to be copied over into the next generation
	 */
	public Prisoner[] fitnessProportionateSelection(Prisoner[] c, int m) {
		popSize = c.length;
		Prisoner Selected[] = new Prisoner[popSize];

		// When elitism is enabled
		if(m != 0) {
			Arrays.sort(c);
			Prisoner Elite[];
			Elite = Arrays.copyOfRange(c, 0, m);
			for(int i = 0; i < m; i++) {
				Selected[i] = Elite[i];
			}
		}

        /* re-evaluate the fitnesses of all individuals using the sigma scaling
           method with the parameter given in the book */

        // calculating the mean fitness
        double totalFitness = 0;
		for(int i = 0; i < popSize; i++) {
		    totalFitness += c[i].getScore();
        }
        double meanFitness = totalFitness / popSize;

		// calculate the variance and standard deviation
        double sumFitnessDiff = 0;
        for(int i = 0; i < popSize; i++) {
            sumFitnessDiff += Math.pow(c[i].getScore() - meanFitness, 2);
        }
        double variance = sumFitnessDiff / popSize;
        double std_deviation = Math.sqrt(variance);

        // set fitness to new scaled fitness
        for(int i = 0; i < popSize; i++) {
            double newFitness = (c[i].getScore() - meanFitness) /
                    (2 * std_deviation);
            c[i].setScore(newFitness);
        }

        /* fitness proportional selection with stochastic universal sampling */
		// calculating the distance between pointers
        double newTotalFitness = 0;
        for(int i = 0; i < popSize; i++) {
            newTotalFitness += c[i].getScore();
        }
        int numOffsprings = popSize - m;
        double pointerDistance = newTotalFitness / numOffsprings;

        // setting a random starting point
        Random r = new Random();
		double startPos = pointerDistance * r.nextDouble();

		// generate N equally spaced pointers
		Prisoner[] individuals = new Prisoner[numOffsprings];
		int index = 0;
		double sum = c[index].getScore();

		for(int i = 0; i < numOffsprings; i++) {
			double curPointer = startPos + i * pointerDistance;
			if(sum >= curPointer) {
				individuals[i] = c[index];
			} else {
				for(++index; index < popSize; index++) {
					sum += c[index].getScore();
					if(sum >= curPointer) {
						individuals[i] = c[index];
						break;
					}
				}
			}
		}

		for(int i = m; i < popSize; i++) Selected[i] = individuals[i-m];
		return Selected;

//		for(int i = 0; i < numOffsprings; i++) {
//			double curPointer = startPos + i * pointerDistance;
//			pointers[i] = curPointer;
//		}
//
//		Prisoner RWSSelected[] = RouletteWheelSelection(c, pointers);
//		for(int i = m; i < popSize; i++) Selected[i] = RWSSelected[i-m];
//		return Selected;
	}

    /**
     * A helper function used for FPS with SUS
     * @param pop
     * @param pointers
     * @return
     */
	public Prisoner[] RouletteWheelSelection(Prisoner[] pop, double[] pointers) {
		int pop_size = pop.length;
		Prisoner Keep[] = new Prisoner[pop_size];
		for (int curPos = 0; curPos < pointers.length; curPos++) {
			double p = pointers[curPos];
			int curPrisoner = 0;
			double fitSum = 0;
			while(fitSum < p) {
				fitSum = 0;
				for (int i = 0; i < curPrisoner; i++) {
					fitSum += pop[curPrisoner].getScore();
				}
				curPrisoner++;
			}
			Keep[curPos] = pop[curPrisoner];
		}
	    return Keep;
    }

    /**
     * Enables tournament selection
	 * @param c - initial population
	 * @param k - number of unique individuals to be chosen from each tournament
     */
    public Prisoner[] tournamentSelection(Prisoner[] c, int k){
    	curPopulation = c;
    	popSize = curPopulation.length;
		Prisoner keep[] = new Prisoner[popSize];
		Prisoner round[] = new Prisoner[k];

		// hosting popSize rounds of tournaments
		for(int i = 0; i < popSize; i++) {
			// for each tournament, choose k unique individuals
			for(int j = 0; j < k ; j++) {
				Random rand = new Random();
				int value = rand.nextInt(popSize);
				while (Arrays.asList(round).contains(curPopulation[value]) && value < popSize-1) {
					value++;
				}
				round[j] = curPopulation[value];
			}
			// after picking k unique individuals, find the individual with the highest fitness
			double maxFit = 0;
			int indexBest = 0;
			for (int m = 0; m < k; m++) {
				if (round[m].getScore() > maxFit) {
					maxFit = round[m].getScore();
					indexBest = m;
				}
			}
			// adding the individual with the highest fitness to the list of individuals to keep
			keep[i] = round[indexBest];
		}
		return keep;
    }

    /**
     *Responsible for updating the graphical representation
     */
    public void paintComponent(Graphics g) {
        super.paintComponent(g); //paint background

		//Get limits of viewable area
      	Insets insets = getInsets();
      	int x0 = insets.left;
      	int y0 = insets.top;
		int currentWidth = getWidth() - insets.left - insets.right;
		int currentHeight = getHeight() - insets.top - insets.bottom;
	    
		//Display a series of rectangles, representing the players
		for(int i = 0; i < popSize; i++) {
	    	g.setColor(curPopulation[i].getColor());	
	    	g.fillRect((x0*2)+((currentWidth/popSize)*(i)),(currentHeight/4)+y0,(currentWidth/popSize),currentHeight/2);
	    }
    }
	
    /**
     *Set the {@link  ie.errity.pd.Rules Rules}
     *@param r new {@link  ie.errity.pd.Rules Rules}
     */
    public void setRules(Rules r) {
		mutateP = r.getMutateP();
		crossP = r.getCrossP();
		selection = r.getSelection();
		selParam = r.getSelectionParameter();
		seed = r.getSeed();
		if (seed > -1) rand.setSeed(seed);
	}

	/**
	 *Reset the breeder
	 */
	public void clear() {
		popSize = 0;
		repaint();
	}
}
	
