/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package org.ualg.pt.pmfcore;

import java.util.ArrayList;
import java.util.Collections;
import java.util.Comparator;
import java.util.HashMap;
import java.util.HashSet;
import java.util.Iterator;
import java.util.LinkedHashMap;
import java.util.LinkedList;
import java.util.List;
import java.util.ListIterator;
import java.util.Map;
import java.util.Set;


/**
 *
 * @author Eman
 */
public class ProteinMatcher {

    double [] proteinFreq;
    double[][] frequencyTable;
    double [][] frequencyTable_RowNormalized;
    int maxProteinMass;
    int maxPeptidesMass;
    final int PROTIEN_INTERVAL = 10000;
    final int PEPTIDE_INTERVAL = 100;

    public ProteinMatcher(int maxProteinMass, int maxPeptidesMass) {
        this.maxProteinMass = maxProteinMass;
        this.maxPeptidesMass = maxPeptidesMass;
        frequencyTable = new double[(maxPeptidesMass / PEPTIDE_INTERVAL) + 1][(maxProteinMass / PROTIEN_INTERVAL) + 1];
        proteinFreq=new double[(maxProteinMass/PROTIEN_INTERVAL)+1];
    }

    /**
     *Get frequency table for a protein
     * @return array of frequency of protein masses values.
     */
    public double[] getProteinFreq() {
        return proteinFreq;
    }

    /**
     * Set frequency table for a protein
     * @param proteinFreq
     */
    public void setProteinFreq(double[] proteinFreq) {
        this.proteinFreq = proteinFreq;
    }
    
    

    /**
     * Get the frequency table of database proteins
     * @return frequency table
     */
    public double[][] getFrequencyTable() {
        return frequencyTable;
    }

    /**
     * Set the frequency table
     * @param frequencyTable
     */
    public void setFrequencyTable(double[][] frequencyTable) {
        this.frequencyTable = frequencyTable;
    }

    /**
     * Get the maximum molecular weight of protein database
     * @return maximum molecular weight
     */
    public int getMaxProteinMass() {
        return maxProteinMass;
    }

    /**
     * Set the maximum molecular weight of protein database
     * @param maxProteinMass
     */
    public void setMaxProteinMass(int maxProteinMass) {
        this.maxProteinMass = maxProteinMass;
    }

    /**
     * Get the maximum  peptide molecular weight in protein database
     * @return maximum  peptide molecular weight
     */
    public int getMaxPeptidesMass() {
        return maxPeptidesMass;
    }

    /**
     * Set the maximum  peptide molecular weight in protein database
     * @param maxPeptidesMass
     */
    public void setMaxPeptidesMass(int maxPeptidesMass) {
        this.maxPeptidesMass = maxPeptidesMass;
    }

   
    /**
     * ***************************************************************************************************************
     * Build the frequency table in which cell contains the number of occurrence
     * within specific peptides mass range in a certain protein mass.
     *
     * @param proteins
     * @param maxProteinMass
     * @param maxPeptidesMass
     */
    public void addFrequencies(Protein protein) {
        if (frequencyTable == null) {
            return;
        }
      
        int width = frequencyTable.length;
        int height = frequencyTable[0].length;
        Iterator massIterator = protein.getMassesList().iterator(); // get masses list of the current protein      
        while (massIterator.hasNext()) {
            // Get current mass then compare it to all other masses in the other proteins
            double currentMass = (Double) (massIterator.next());
            frequencyTable[(int) currentMass / PEPTIDE_INTERVAL][(int) (protein.getProteinMass() / PROTIEN_INTERVAL)]++;
        }
    }
    
    
    /**
     * Add frequencies to table, in which each column (represents protein mass) is divided into 10000 Dalton, 
     * and each row  (represents peptide mass)is divided into 100 Dalton 
     * @param Piptidemasses
     * @param proteinMass
     */
    public void addFrequencies(Set<Double> Piptidemasses, double proteinMass){
         
        for(Double mass:Piptidemasses){
           frequencyTable[(int) (mass / PEPTIDE_INTERVAL)][(int) (proteinMass / PROTIEN_INTERVAL)]++;
        }
    }
    
    
    /**
     * ***************************************************************************************************************
     * Print frequency table each column is 10 KDa interval and each row is 100
     * Dalton interval
     */
    public void printFrequenceyTable() {
        for (int i = 0; i < frequencyTable.length; i++) {
            for (int j = 0; j < frequencyTable[0].length; j++) {
                System.out.print(frequencyTable[i][j] + "\t\t");
            }
            System.out.println();
        }
    }
    
    
    /**
     * *******************************************************************************************************************
     * get the largest value in each column
     *
     * @return arrayList of integer
     */
    public ArrayList<Double> largestValuesInColumns() {

        ArrayList<Double> maximumColumnsValues = new ArrayList<Double>();
        double maxValueInColumn;
        for (int i = 0; i < frequencyTable[0].length; i++) {
            maxValueInColumn = Double.MIN_VALUE;

            for (int j = 0; j < frequencyTable.length; j++) {
                if (frequencyTable[j][i] > maxValueInColumn) {
                    maxValueInColumn = frequencyTable[j][i];
                }
            }
            maximumColumnsValues.add(maxValueInColumn);
        }
        return maximumColumnsValues;
    }

    /**
     * **************************************************************************************************************
     * Normalize every cell value in the table to the largest value in each
     * column.
     */
    public void normalizeToLargestValue() {

        ArrayList<Double> largestColValues = largestValuesInColumns();
        // DecimalFormat deFormat=new DecimalFormat("#.####");
        for (int i = 0; i < frequencyTable[0].length; i++) {
            for (int j = 0; j < frequencyTable.length; j++) {
                if (largestColValues.get(i) == 0.0) {
                    continue;
                }
                frequencyTable[j][i] = ((frequencyTable[j][i] / largestColValues.get(i)));
            }
        }
    }
 
    /**
     * ***************************************************************************************************************
     * Compare experimental masses with each protein masses in the database and
     * for each match return the frequency from frequency table
     *
     * @param exProtein
     * @param thProtein
     * @return list of matched proteins masses
     */
    public void doMatching(Protein exProtein, ArrayList<Protein> thProtein) {
        int highestScoresNumber=0;
        double result = 0;
        double currentThMass;
        double currentExMass;
        double threshold = 1;
        ArrayList<Double> peptideFreq = new ArrayList<Double>();
        Map<Protein, Double> occuranceMap = new HashMap<Protein, Double>();
        Map<Protein, Double> sortedOccuranceMap = new HashMap<Protein, Double>();
        Iterator thIterator = thProtein.iterator();
        while (thIterator.hasNext()) {
            // Get the current theoretical protein from the database
            Protein thP = (Protein) thIterator.next();
            Iterator massExIterator = exProtein.getMassesList().iterator();// get masses list of the experimental protein
            ArrayList<Double> matchedMasses = new ArrayList<Double>();

            while (massExIterator.hasNext()) {
                currentExMass = (Double) massExIterator.next(); // Get current mass then compare it to all other masses in the other proteins
                Iterator massThIterator = thP.getMassesList().iterator(); // Get the masses list of this protein
                while (massThIterator.hasNext()) {
                    currentThMass = (Double) massThIterator.next();
                    // Do the comparision
                    result = Math.abs(currentThMass - currentExMass);
                    if (result <= threshold) {
                        matchedMasses.add(currentThMass);
                    }
                }
            }
            
            if (!matchedMasses.isEmpty()) {
                // Sort the map to get the highest 10 scores
                sortedOccuranceMap = sortByValue(occuranceMap);
            }
        }
        
        //Reverse the map and print the top scores. 
        ListIterator<Map.Entry<Protein, Double>> iter =
                new ArrayList(sortedOccuranceMap.entrySet()).listIterator(sortedOccuranceMap.size());
        
        while (iter.hasPrevious()) {
            Map.Entry<Protein, Double> ent = iter.previous();
            Protein p = (Protein) ent.getKey();
            double scoreValue = (Double) ent.getValue();
            System.out.println(" Its Score:" + scoreValue + "  Protein mass: " + p.getProteinMass() + "  Protein masses list: " + p.getMassesList());
            highestScoresNumber++;
        }
    }
    
    
    /**
     * ***************************************************************************************************************
     * Compare experimental masses with each protein masses in the database and
     * for each match return the frequency from frequency table
     *
     * @param exProtein Protein Sample: Experimental Masses
     * @param thProtein Database Protein: Theoretical Masses
     * @return list of matched proteins masses
     */
    public double doMatchingOneByOne(Protein exProtein, Protein thP, double threshold) {
        
        double result = 0;
        double score = 0;
        double currentThMass;
        double currentExMass;
        thP.numberOfMatched=0;
        ArrayList<Double> peptideFreq = new ArrayList<Double>();
        Map<Protein, Double> occuranceMap = new HashMap<Protein, Double>();
        Map<Protein, Double> sortedOccuranceMap = new HashMap<Protein, Double>();
            Iterator massExIterator = exProtein.getMassesList().iterator();// get masses list of the experimental protein
            Set<Double> matchedMasses = new HashSet<Double>();
   
            while (massExIterator.hasNext()) {
                currentExMass = (Double) massExIterator.next(); // Get current mass then compare it to all other masses in the other proteins
                if(thP.getMassesList().isEmpty()) {
                    System.out.println("theoretical Protein  has no masses: " + thP.ProteinID);
                    continue;
                }
                for(Peptide thpeptide: thP.getPeptides()){
                    currentThMass =  thpeptide.getMass();
                    
                    // Do the comparision
                  // if(thpeptide.getPeptideLength()<3) continue; // ignore small piptides
                   result = Math.abs(currentThMass - currentExMass);
                    
                    if (result <= threshold) { // matched
                        matchedMasses.add(currentThMass);
                        thP.addMatchedPeptide(thpeptide);
                        thP.numberOfMatched++;
                    }
                }
            }
                      
            if (!matchedMasses.isEmpty()) {
                //Here Call method to score theoretical protein and pass the list of matched masses
                 score = getPeptidesFrequenceyForProtein(thP, matchedMasses);
              
            }          
             score= score*thP.numberOfMatched*((double)thP.numberOfMatched/(double)exProtein.getMassesList().size())*((double)thP.getMatchedCoverage());      
              return score;
    }
    
    /**
     * Assign the scores for the  proteins Database and sort them from the highest to the lowest.
     * @param ExProtein
     * @param thProteins
     */
    public void proteinsScores(Protein ExProtein,ArrayList<Protein> thProteins, double dataBaseThreshold){
        Map<Protein, Double> sortedOccuranceMap = new HashMap<Protein, Double>();
        Iterator dbProteinsIt=thProteins.iterator();
        Map<Protein, Double> occuranceMap = new HashMap();
        int highestScoresNumber=0;
        while(dbProteinsIt.hasNext()){
            Protein oneProtein=(Protein) dbProteinsIt.next();
            double score=doMatchingOneByOne(ExProtein, oneProtein,dataBaseThreshold);
            occuranceMap.put(oneProtein, score);     
        }
        sortedOccuranceMap = sortByValue(occuranceMap);
        ListIterator<Map.Entry<Protein, Double>> iter =
                new ArrayList(sortedOccuranceMap.entrySet()).listIterator(sortedOccuranceMap.size());
        
        while (iter.hasPrevious()) {
            Map.Entry<Protein, Double> ent = iter.previous();
            Protein p = (Protein) ent.getKey();
            double scoreValue = (Double) ent.getValue();
            System.out.println("\nProtein ID: "+p.getProteinID()+" Its Score:" + scoreValue + "  Protein mass: " + p.getProteinMass());
            highestScoresNumber++;
            //to test just print the top 3 scores
            if(highestScoresNumber>5) break;
        }       
    }

    /**
     * ***************************************************************************************************************
     * Get the frequency of peptides for each protein
     *
     * @param matchedPeptides
     * @param thprotein
     */
    public double getPeptidesFrequenceyForProtein(Protein thprotein, Set<Double> matchedMasses) {
        double occurance = 1.0;
        double score = 0.0;
        double proteinMass = thprotein.getProteinMass();
        int protien_mass_index_in_frequency_table = (int) (proteinMass / PROTIEN_INTERVAL);
        Iterator it = matchedMasses.iterator();
        Set<Integer> validMasses= new HashSet<Integer>();
        while (it.hasNext()) {
            double mass = (Double) it.next();
            int mass_index_in_frequency_table = (int) (mass / PEPTIDE_INTERVAL);
            if(occurance==0) {
                score=0.0;
            }
            occurance *= frequencyTable[mass_index_in_frequency_table][protien_mass_index_in_frequency_table];
        }
      
        occurance=occurance/(thprotein.numberOfMatched+1);
        score =  50000/(occurance*proteinMass);
        return score;
    }

    
    /**
     * Sort the map that contains the database proteins scores
     *
     * @param map
     * @return
     */
    public Map sortByValue(Map map) {
        List list = new LinkedList(map.entrySet());
        Collections.sort(list, new Comparator() {
            @Override
            public int compare(Object o1, Object o2) {
                return ((Comparable) ((Map.Entry) (o1)).getValue())
                        .compareTo(((Map.Entry) (o2)).getValue());
            }
        });

        Map result = new LinkedHashMap();
        for (Iterator it = list.iterator(); it.hasNext();) {
            Map.Entry entry = (Map.Entry) it.next();
            result.put(entry.getKey(), entry.getValue());
        }
        return result;
    } 
        
    /*
     * Get the number of proteins in each interval.
     */
    public void addFrequenciesforProteins(Protein protein) {

        if (proteinFreq == null) {
            return;
        }

        double currentProteinMass = protein.getProteinMass();
        proteinFreq[(int) currentProteinMass / PROTIEN_INTERVAL]++;       
    }

    /**
     * Print the frequency of masses values for protein
     */
    public void printProteinFrequences() {

        for (int i = 0; i < proteinFreq.length; i++) {
            System.out.print(proteinFreq[i] + "\t");
        }
        System.out.println();
    }
   
}
