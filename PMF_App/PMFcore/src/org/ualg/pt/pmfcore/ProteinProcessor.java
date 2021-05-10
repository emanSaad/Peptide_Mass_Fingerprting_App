/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package org.ualg.pt.pmfcore;

import java.io.File;
import java.io.FileNotFoundException;
import java.io.FileOutputStream;
import java.io.IOException;
import java.io.PrintWriter;
import java.util.ArrayList;
import java.util.Collections;
import java.util.HashMap;
import java.util.Iterator;
import java.util.LinkedList;
import java.util.Map;
import java.util.Set;
import java.util.logging.Level;
import java.util.logging.Logger;
import org.biojava.bio.BioException;
import org.biojava.bio.proteomics.Protease;
import org.biojava.bio.seq.ProteinTools;
import org.biojava.bio.symbol.Alphabet;
import org.biojava.bio.symbol.AlphabetManager;
import org.biojava.bio.symbol.Edit;
import org.biojava.bio.symbol.FiniteAlphabet;
import org.biojava.bio.symbol.IllegalAlphabetException;
import org.biojava.bio.symbol.IllegalSymbolException;
import org.biojava.bio.symbol.Symbol;
import org.biojava.utils.ChangeVetoException;
import org.biojavax.bio.seq.RichSequence;



/**
 *
 * @author Eman
 */
public class ProteinProcessor {

    public void doAllThProteinProcess(Protein thprotein) {
    }

    /**
     * *************************************************************************************************
     * Read all the masses files from specific directory
     *
     * @param dir
     * @return
     * @throws FileNotFoundException
     * @throws IOException
     */
    public ArrayList<Protein> readMassesFilesFromDir(File dir) throws FileNotFoundException, IOException {

        File fList[] = dir.listFiles();

        ArrayList<Protein> proteins = new ArrayList<Protein>();
        for (int i = 0; i < fList.length; i++) {
            if (fList[i].getName().equals(".DS_Store") || !(fList[i].isFile())) {
                continue;
            }

            ArrayList<Double> massesValuse = DataReader.readMasses(fList[i]);

            // create a new experemental protien
            Protein protein = new Protein();
            protein.setMassesList(massesValuse);
            protein.ProteinID = fList[i].getName();
            protein.setFile(fList[i]);
            protein.setSequence(null);
            protein.setType(Protein.ProteinType.EXPERIMENTAL);

            // add it to the list
            proteins.add(protein);
        }
        return proteins;

    }
    
   
    /**
     * Read MS masses from one file
     * @param fileName
     * @return
     * @throws FileNotFoundException
     * @throws IOException
     */
    public Protein readMassesFile(File fileName) throws FileNotFoundException, IOException {
        
            ArrayList<Double> massesValuse = DataReader.readMasses(fileName);
            // create a new experemental protien
            Protein protein = new Protein();
            protein.setMassesList(massesValuse);
            protein.ProteinID = fileName.getName();
            protein.setFile(fileName);
            protein.setSequence(null);
            protein.setType(Protein.ProteinType.EXPERIMENTAL);

        return protein;

    }
    
   
    
    
    /**
     * *************************************************************************************************
     * Read sequences and create a list of the ٍprotein sequences
     *
     * @param dir
     * @return a list of theoretical proteins
     * @throws FileNotFoundException
     * @throws IOException
     */
    public ArrayList<Protein> readSequencesFilesFromDir(File dir) throws FileNotFoundException, IOException {

        File fList[] = dir.listFiles();
        ArrayList<Protein> proteins = new ArrayList<Protein>();
        for (int i = 0; i < fList.length; i++) {

            if (fList[i].getName().equals(".DS_Store")) {
                continue;
            }
            System.out.println("Loading Sequence file:" + fList[i].getAbsolutePath());
            ArrayList<RichSequence> sequences = DataReader.readSequence(fList[i]);
            Iterator it = sequences.iterator();
            while (it.hasNext()) {
                // in case of more protein sequences in the same file we create a protein for each sequence.
                // create a new theortical protien
                RichSequence sequence = (RichSequence) it.next();
                Protein protein = new Protein();
                protein.setMassesList(null);
                protein.ProteinID = sequence.getAccession();
                protein.setFile(fList[i]);
                protein.setSequence(sequence);
                protein.setType(Protein.ProteinType.THEORETICAL);
                // add it to the list
                proteins.add(protein);
            }
        }
        return proteins;
    }
    

    /**
     * *************************************************************************************************
     * Digest the protein and calculate the masses for protein by calling method
     * in Digester class
     *
     * @param protiens
     * @throws BioException
     */
    public void digestProteins(ArrayList<Protein> protiens, boolean fixedModification, int maxMissedCleavages) throws BioException {
        Iterator it = protiens.iterator();
        while (it.hasNext()) {
            Protein protein = (Protein) it.next();
            
            protein.setMassesList(DigesterAndMassCalculater.proteinDigestionAndMassesCalculation1(protein, Protease.TRYPSIN,fixedModification, maxMissedCleavages));
            int i=0;
        }
    }

    /**
     * *************************************************************************************************
     * Digest one protein at time and calculate the its peptides masses.
     *
     * @param protein
     */
    public void digestOneProtein(Protein protein,boolean fixedModification, int maxMissedCleavages)  {
        try {
            
            protein.setMassesList(DigesterAndMassCalculater.proteinDigestionAndMassesCalculation1(protein, Protease.TRYPSIN,fixedModification, maxMissedCleavages));
        } catch (BioException ex) {
            Logger.getLogger(ProteinProcessor.class.getName()).log(Level.SEVERE, null, ex);
        } catch (ChangeVetoException ex) {
            Logger.getLogger(ProteinProcessor.class.getName()).log(Level.SEVERE, null, ex);
        }
    }

    /**
     * *************************************************************************************************
     * Print a list of proteins
     *
     * @param protiens
     */
    public void printProtiens(ArrayList<Protein> protiens) {
        Iterator it = protiens.iterator();
        //System.out.println(protiens);

        while (it.hasNext()) {
            Protein protein = (Protein) it.next();
            System.out.println(protein);
            System.out.println("\n");
        }
    }

    /**
     * *************************************************************************************************
     * Print one protein at time because we will process them (proteins)
     * separately.
     *
     * @param protein
     */
    public void printOneProtein(Protein protein) {
        System.out.println(protein);
        System.out.println("\n");
    }

    

    /**
     * ***********************************************************************************************
     * Calculate the masses for input proteins
     *
     * @param proteins
     * @return an array of protein masses
     */
    public ArrayList<Double> calcProteinMass(ArrayList<Protein> proteins) {

        ArrayList<Double> massesProtein = new ArrayList<Double>();
        double proteinMass;
        Iterator it = proteins.iterator();

        while (it.hasNext()) {
            Protein protein = (Protein) it.next();
            RichSequence sequence = protein.getSequence();

            try {
                proteinMass = DigesterAndMassCalculater.getProteinMass(sequence);
                protein.setProteinMass(proteinMass);
                massesProtein.add(proteinMass);
            } catch (IllegalSymbolException ex) {
                Logger.getLogger(ProteinProcessor.class.getName()).log(Level.SEVERE, null, ex);
            }
        }
        return massesProtein;
    }
    
    /**
     * Calculate the molecular weight for one protein
     * @param protein
     * @return
     * @throws IllegalSymbolException
     */
    public double calcProteinMass(Protein protein) throws IllegalSymbolException {
       
        double proteinMass=0.0;
            RichSequence sequence = protein.getSequence();

                proteinMass = DigesterAndMassCalculater.getProteinMass(sequence);
                protein.setProteinMass(proteinMass);
            
        //System.out.println("Intact protein mass:\n" + massesProtein );
        return proteinMass; 
}
    

    /**
     * *************************************************************************************************
     * get the maximum protein mass(molecular weight)
     *
     * @param proteinMasses
     * @return the largest protein mass
     */
    public double getMaximumProteinMass(ArrayList<Double> proteinMasses) {

        double maxProteinMass = Collections.max(proteinMasses);

        return maxProteinMass;
    }

    /**
     * *************************************************************************************************
     * give the minimum value of protein masses
     *
     * @param proteinMasses
     * @return lowest protein mass
     */
    public double getMinimumProteinMass(ArrayList<Double> proteinMasses) {

        double minProteinMass = Collections.min(proteinMasses);
        return minProteinMass;
    }

    /**
     * *************************************************************************************************
     * Get the maximum peptide mass
     *
     * @param proteins
     * @return the largest peptide mass
     */
    public double getMaximumPeptideMass(ArrayList<Protein> proteins) {
        double maxPeptideMass=0.0;
        Iterator iterator = proteins.iterator();
        
        while (iterator.hasNext()) {
            Protein protein = (Protein) iterator.next();
            double currentMax = Collections.max(protein.getMassesList());
            if (maxPeptideMass < currentMax) {
                maxPeptideMass = currentMax;
            }
        }

        return maxPeptideMass;
    }

    /**
     * *************************************************************************************************
     * Remove the contaminants from the experimental masses with specified
     * threshold by digest the known Protein contaminants and calculate the
     * masses for them. Then, compare the Resulting masses to experimental
     * masses
     *
     * @param contProtein
     * @param exProtein
     * @param threshold
     * @return ArrayList of valid masses
     */
    public ArrayList<Double> massesContaminantsRemoval(ArrayList<Protein> exProtein, ArrayList<Protein> contProtein, double threshold) {
        ArrayList<Double> validMasses = new ArrayList<Double>();
        ArrayList<Double> contMasses = new ArrayList<Double>();

        double result = 0;
        
        Iterator exIterator = exProtein.iterator();
        while (exIterator.hasNext()) {
            Protein exP = (Protein) exIterator.next(); // get the current protein

            Iterator massExIterator = exP.getMassesList().iterator(); // get masses list of the current protein
            while (massExIterator.hasNext()) {
                // Get current mass then compare it to all other masses in the other proteins
                double currentExMass = (Double) massExIterator.next();

                Iterator contIterator = contProtein.iterator();
                while (contIterator.hasNext()) {
                    Protein contaminantP = (Protein) contIterator.next();

                    Iterator massContIterator = contaminantP.getMassesList().iterator();
                    while (massContIterator.hasNext()) {
                        double currentcontMass = (Double) massContIterator.next();
                        // Do the comparision
                        result = Math.abs(currentcontMass - currentExMass);
                        if (result <= threshold) {
                            // add the experimental masses that are appear to be contaminants to contaminants list
                            contMasses.add(currentExMass);                                                 
                        }
                    }
                }
                if (result > threshold) {
                    System.out.println(currentExMass);
                    validMasses.add(currentExMass);
                }

            }
        }

        System.out.println("Valid experimental masses: " + validMasses);
        System.out.println("Contaminants masses: " + contMasses);
        return validMasses;
    }
    
    

    /**
     * **************************************************************************************************
     * Write Valid experimental masses to file
     *
     * @param validMasses
     * @param validMassesfile
     * @throws FileNotFoundException
     */
    public void writeValidMassesToFile(ArrayList<Double> validMasses, String validMassesfile) throws FileNotFoundException {
        PrintWriter pw = new PrintWriter(new FileOutputStream(validMassesfile));
        Iterator it = validMasses.iterator();
        while (it.hasNext()) {
            double validMass = (Double) it.next();
            pw.println(validMass);
        }
    }

    
    
    /**
     * *************************************************************************************************
     * get the protein symbols and return them as array list
     *
     * @return protein symbols list
     */
    public ArrayList<Symbol> getProteinSymbol() {
        ArrayList<Symbol> proteinSymbols = new ArrayList<Symbol>();
        FiniteAlphabet protein = ProteinTools.getAlphabet();
        Iterator proteinSym = protein.iterator();
        while (proteinSym.hasNext()) {
            Symbol s = (Symbol) proteinSym.next();
            if (s.getName().equals("SEC") || s.getName().equals("PYL")) {
               continue;
            }
            proteinSymbols.add(s);
            //System.out.println("List of protein symbols: " + s.getName());
        }
        return proteinSymbols;
    }

    /**
     * *************************************************************************************************
     * Get the X in sequence and replace it with all available amino acids
     * codes, so one X in sequence it needs to 22 combination when the sequence
     * has 2 X, it needs to 22^2. so we need to (22)^number of Xs.
     *
     * @param thProtein
     * @param symList
     */
    public void translateXtoProteinSymbol(ArrayList<Protein> thProteins) throws IllegalSymbolException, IllegalAlphabetException, BioException {

        ArrayList<Protein> newList = new ArrayList<Protein>();  // a new arralist where 
        ArrayList<Protein> todelete = new ArrayList<Protein>();
        Iterator itp = thProteins.iterator();
        ArrayList<Symbol> listofSymbols = getProteinSymbol();
        System.out.println("Number of symbols:" + listofSymbols.size());
        ProteinChecker proteinChecker = new ProteinChecker();
        while (itp.hasNext()) {
            Protein thProtein = (Protein) itp.next();

            Map<Integer, Symbol> unknownSymbols = proteinChecker.getUnknowSymbols(thProtein, true);
            if (unknownSymbols.isEmpty()) {
                continue; // no unknown symbol was found.
            }
            int count = unknownSymbols.size(); // number of unknown symbols
            if (count > 4) {
                System.out.println("Too many unknown symbols. Quiting!");
                System.exit(0);
            }

            /*       Compensation of the unknown symbols
             *   The totall number of new proteins will be = proteinSymbols^count
             *   For each compensation symbol create a list of proteins, each one of them (protein) has this symbol.
             *   add this list to the protein list and go ahead.
             */

            Iterator unknownSymbolsIterator = unknownSymbols.entrySet().iterator();
            if (unknownSymbolsIterator.hasNext()) { // get only the first unknow symbol the propogate 
                // get current entry
                Map.Entry entry = (Map.Entry<Integer, Symbol>) unknownSymbolsIterator.next();
                int position = (Integer) entry.getKey();
                Symbol symP = (Symbol) entry.getValue();


                // get it's name if exists, otherwise show an error.
                String name = symP.getName();
                if (name == null) { // show error and return
                    System.out.println("[INF] In protein:" + thProtein.getProteinID() + "\t (unknown symbol)");
                    continue;
                }

                System.out.println("Size of Protein:" + thProtein.getSequence().toList().size() * 2);
                LinkedList<Protein> queue = new LinkedList<Protein>();

                // add all possible symbols to the queue.
                for (int i = 0; i < listofSymbols.size(); i++) {
                    Edit edRest = new Edit(position, thProtein.getSequence().getAlphabet(), listofSymbols.get(i));
                    Protein newProtein = new Protein(thProtein);
                    newProtein.getSequence().edit(edRest);
                    queue.add(newProtein);
                }

                // now remove the protein
                todelete.add(thProtein);
                // process the queue
                while (!queue.isEmpty()) {
                    Protein current = queue.pollFirst(); // we need to test it first.
                    // if the current protein has no more X symbols, then just remove it and add it to the new list
                    Map<Integer, Symbol> unknownsymbols = proteinChecker.getUnknowSymbols(current, false);
                    if (unknownsymbols.isEmpty()) { // well it is clear add it to the newlist

                        newList.add(current);
                        continue;
                    }                   
                    // get the first entry in the map
                    // we are sure we have at least one item in the map, so we need not to check with .hasnext function
                    // add new proteins to the queue

                    Map.Entry ent = (Map.Entry<Integer, Symbol>) unknownsymbols.entrySet().iterator().next();
                    int pos = (Integer) ent.getKey();
                    for (int i = 0; i < listofSymbols.size(); i++) {
                        Edit edRest = new Edit(pos, current.getSequence().getAlphabet(), listofSymbols.get(i));
                        Protein newProtein = new Protein(current);
                        newProtein.getSequence().edit(edRest);
                        queue.add(newProtein);
                    }
                }
            }

        }
        // add the new proteins to to main list;
        thProteins.addAll(newList);
        // delete the Xed proteins
        thProteins.removeAll(todelete);

    }


    /**
     * ****************************************************************************************************
     * Get an array contains the positions of Xs in the protein sequence
     *
     * @param thProtein
     * @return array
     */
    public ArrayList<Integer> getPositionsOfX(Protein thProtein) {

        ArrayList<Integer> positionsForXs = new ArrayList<Integer>();
        int unknowSymbolPosition = 1;
        RichSequence sequence = thProtein.getSequence();
        Iterator sequenceIt = sequence.iterator();
        while (sequenceIt.hasNext()) {
            Symbol symbol = (Symbol) sequenceIt.next();
            if (symbol.getMatches().toString() == null) {
                positionsForXs.add(unknowSymbolPosition);

            }
            unknowSymbolPosition++;
        }
       
        return positionsForXs;
    }
    
    /**
     * ****************************************************************************************************
     * Get an array contains the positions of Us in the protein sequence
     *
     * @param thProtein
     * @return array
     */
    public ArrayList<Integer> getPositionsOfU(Protein thProtein) {

        ArrayList<Integer> positionsForUs = new ArrayList<Integer>();
        int uSymbolPosition = 1;
        RichSequence sequence = thProtein.getSequence();
        Iterator sequenceIt = sequence.iterator();
        while (sequenceIt.hasNext()) {
            Symbol symbol = (Symbol) sequenceIt.next();
            if (symbol.getName().equals("SEC")) {
                positionsForUs.add(uSymbolPosition);

            }
            uSymbolPosition++;
        }
        //System.out.println("positionOf Xs ar : " + unknowSymbolPosition);
        return positionsForUs;
    }
    
    
    /**
     * ****************************************************************************************************
     * Get an array contains the positions of Os in the protein sequence
     *
     * @param thProtein
     * @return array
     */
    public ArrayList<Integer> getPositionsOfO(Protein thProtein) {

        ArrayList<Integer> positionsForOs = new ArrayList<Integer>();
        int oSymbolPosition = 1;
        RichSequence sequence = thProtein.getSequence();
        Iterator sequenceIt = sequence.iterator();
        while (sequenceIt.hasNext()) {
            Symbol symbol = (Symbol) sequenceIt.next();
            if (symbol.getName().equals("PYL")) {
                positionsForOs.add(oSymbolPosition);

            }
            oSymbolPosition++;
        }
       // System.out.println("positionOf Os ar : " + oSymbolPosition);
        return positionsForOs;
    }
    
    

    /**
     * ****************************************************************************************************
     * Get an array of the positions of unknown symbols in the protein sequence
     *
     * @param thProtein
     * @return array
     */
    public Map<String, ArrayList<Integer>> getPositionsOfUnknownSymbols(Protein thProtein) throws IllegalSymbolException {

        Map<String, ArrayList<Integer>> positionsForUnknowns = new HashMap<String, ArrayList<Integer>>();
        FiniteAlphabet protein = ProteinTools.getAlphabet();
        ArrayList<Integer> positionsofX = new ArrayList<Integer>();
        ArrayList<Integer> positionsofB = new ArrayList<Integer>();
        ArrayList<Integer> positionsofZ = new ArrayList<Integer>();
        ArrayList<Integer> positionsofJ = new ArrayList<Integer>();

        Symbol symbolX;
        Symbol symbolB;
        Symbol symbolZ;
        Symbol symbolJ;
        int Positions = 1;


        RichSequence sequence = thProtein.getSequence();
        Iterator sequenceIt = sequence.iterator();
        Set<Symbol> symbolofX= DigesterAndMassCalculater.getXsymbols();
        symbolX=protein.getAmbiguity(symbolofX);

        Set<Symbol> symbolsOfB = DigesterAndMassCalculater.getBsymbols();
        symbolB = protein.getAmbiguity(symbolsOfB);

        Set<Symbol> symbolsOfZ = DigesterAndMassCalculater.getZsymbols();
        symbolZ = protein.getAmbiguity(symbolsOfZ);

        Set<Symbol> symbolsOfJ = DigesterAndMassCalculater.getJsymbols();
        symbolJ = protein.getAmbiguity(symbolsOfJ);

        while (sequenceIt.hasNext()) {

            Symbol symbol = (Symbol) sequenceIt.next();

            if (symbol.getMatches().contains(symbolX) && symbol != symbolB && symbol != symbolZ && symbol != symbolJ) {
                positionsofX.add(Positions);
                positionsForUnknowns.put(symbol.getName().toString(), positionsofX);

            }

            if (symbol.getMatches().contains(symbolB)&& symbol !=symbolX) {
                positionsofB.add(Positions);
                positionsForUnknowns.put(symbol.getName().toString(), positionsofB);
                
            }

            if (symbol.getMatches().contains(symbolZ) && symbol !=symbolX) {
                positionsofZ.add(Positions);
                positionsForUnknowns.put(symbol.getName().toString(), positionsofZ);
            }

            if (symbol.getMatches().contains(symbolJ) && symbol !=symbolX) {
                positionsofJ.add(Positions);
                positionsForUnknowns.put(symbol.getName().toString(), positionsofJ);
            }
            Positions++;
            
        }
        System.out.println("**************************Unknown Symbols ant Their Positions ********************************");
        Iterator mapIt=positionsForUnknowns.entrySet().iterator();
        while(mapIt.hasNext()){
            Map.Entry mapEntry=(Map.Entry)mapIt.next();
            System.out.println("  Symbol Name: " + mapEntry.getKey()+"  Positions:  " + mapEntry.getValue());
        }        
        return positionsForUnknowns;
    }

    
    /**
     * ************************************************************************************************************
     * Count how many Xs in a sequence because if no X, the sequence is a valid,
     * if just one X in the sequence, we will replace the X with 22 amino acid
     * codes generating with that 22 sequences. if the Xs more than one we will
     * modified the each X when calculating the masses with amino acid average
     * mass 110Da.
     *
     * @param sequence
     * @return the number of Xs
     */
    public int checkHowmanyXs(RichSequence sequence) {
        int countTheXs = 0;
        Iterator sequenceIt = sequence.iterator();
        while (sequenceIt.hasNext()) {
            Symbol oneSymbol = (Symbol) sequenceIt.next();
            if (oneSymbol.equals('X')) {
                countTheXs++;
            }
        }
        return countTheXs;
    }

    /**
     * ****************************************************************************************************************
     * Preprocess the Xs that are more one in one sequence by digesting the
     * sequence and put the X mass to 110Da
     *
     * @param thProteins
     * @throws IllegalSymbolException
     * @throws BioException
     */
    public void digestionAndCalculateMassesWithModification(Protein protein, boolean fixedModification, boolean variableModification, int maxMissedCleavages) throws IllegalSymbolException, BioException {
        RichSequence sequence=protein.getSequence();
        protein.setMassesList(DigesterAndMassCalculater.proteinDigestionAndMassesCalculationWithModification(protein, Protease.TRYPSIN, fixedModification, variableModification, maxMissedCleavages));
    }
    
    
    /**
     * *****************************************************************************************************
     * Replace the unknown symbols in the protein by their corresponding known symbols. Get the positions of all Us and
     * start replace them
     *
     * @param thProtein
     * @param positions
     */
    public RichSequence replacementUsWithOneCode(Protein thProtein,ArrayList<Integer> uPositions ) throws IllegalSymbolException, BioException {
         
        String sequenceAfterReplacement ="";
        Alphabet proteinAlpha = AlphabetManager.alphabetForName("PROTEIN");
        String sequenceString = thProtein.getSequence().seqString();
     
        int currentPosition = 0;
        for (int i = 0; i < uPositions.size(); i++) {
            currentPosition = sequenceString.charAt(uPositions.get(i)-1);
            sequenceAfterReplacement =sequenceString.replace((char) currentPosition, (char) 65);                                     
        }
       
        RichSequence sequenceAfter = RichSequence.Tools.createRichSequence("seqi", sequenceAfterReplacement, proteinAlpha);
        thProtein.setSequence(sequenceAfter);
        return sequenceAfter;
    }
    
   
     /**
     * *****************************************************************************************************
     * Replace the unknown symbols in the protein by their corresponding known symbols. Get the positions of all Us and
     * start replace them
     *
     * @param thProtein
     * @param positions
     */
    public RichSequence replacementOsWithOneCode(Protein thProtein,ArrayList<Integer> oPositions ) throws IllegalSymbolException, BioException {
        
        String sequenceAfterReplacement ="";
        Alphabet proteinAlpha = AlphabetManager.alphabetForName("PROTEIN");
        String sequenceString = thProtein.getSequence().seqString();
        int currentPosition;
        
        for (int i = 0; i < oPositions.size(); i++) {
            currentPosition = sequenceString.charAt(oPositions.get(i)-1);
            sequenceAfterReplacement =sequenceString.replace((char) currentPosition, (char) 65);
          
        }
        RichSequence sequenceAfter = RichSequence.Tools.createRichSequence("seqi", sequenceAfterReplacement, proteinAlpha);
        thProtein.setSequence(sequenceAfter);
        return sequenceAfter;
    }
    
    
}
