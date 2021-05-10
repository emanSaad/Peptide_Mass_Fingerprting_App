/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package org.ualg.pt.pmfcore;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileNotFoundException;
import java.io.FileOutputStream;
import java.io.FileReader;
import java.io.IOException;
import java.io.PrintWriter;
import java.util.ArrayList;
import java.util.Collections;
import java.util.HashMap;
import java.util.Iterator;
import java.util.LinkedHashSet;
import java.util.LinkedList;
import java.util.Map;
import java.util.Random;
import java.util.Set;
import java.util.logging.Level;
import java.util.logging.Logger;
import org.biojava.bio.BioException;
import org.biojava.bio.seq.ProteinTools;
import org.biojava.bio.symbol.Alphabet;
import org.biojava.bio.symbol.AlphabetManager;
import org.biojava.bio.symbol.Edit;
import org.biojava.bio.symbol.IllegalAlphabetException;
import org.biojava.bio.symbol.IllegalSymbolException;
import org.biojava.bio.symbol.Symbol;
import org.biojava.utils.ChangeVetoException;
import org.biojavax.SimpleNamespace;
import org.biojavax.bio.seq.RichSequence;
import org.biojavax.bio.seq.RichSequenceIterator;

/**
 *
 * @author Eman
 */
public class ProteinChecker {

    /**
     * Get the positions of the unknown symbols if found
     *
     * @param thProtein
     * @param verbose true show more information during checking the protein
     * @return a map of symbols and positions to the unknown symbols, empty list
     * if the protein has no unknown symbol
     */
    public Map<Integer, Symbol> getUnknowSymbols(Protein thProtein, boolean verbose) {
        Map<Integer, Symbol> symbols = new HashMap<Integer, Symbol>();
        RichSequence rs = thProtein.getSequence();

        if (rs == null) {
            if (verbose) {
                System.out.println("Protein:" + thProtein.getProteinID() + "\thas no symbols");
            }
            return symbols;
        }
        Iterator symbolsIt = rs.iterator();
        int symbolIndex = 1;
        while (symbolsIt.hasNext()) {
            //List l=Collections.nCopies(22, rs);
            Symbol symP = (Symbol) symbolsIt.next();
            if (symP.getMatches().toString() == null || symP==ProteinTools.o() || symP== ProteinTools.u()) { // if it has no matching, it is a unknown symbol "X"
                if (verbose) {
                    System.out.println("[INFO] An unknown symbol was detected in Protein:" + thProtein.getProteinID() + "\t At position:" + symbolIndex);
                }
                symbols.put(symbolIndex, symP);
            }
            symbolIndex++;
        }
        return symbols;
    }
    
    
    public static Map<Symbol, LinkedList<Integer>> getUnknowSymbolsMap(Protein thProtein, boolean verbose) {
        Map<Symbol, LinkedList<Integer>> symbols = new HashMap<Symbol, LinkedList<Integer>>();
        RichSequence rs = thProtein.getSequence();


        if (rs == null) {
            if (verbose) {
                System.out.println("Protein:" + thProtein.getProteinID() + "\thas no symbols");
            }
            return symbols;
        }
        Iterator symbolsIt = rs.iterator();
        int symbolIndex = 1;
        while (symbolsIt.hasNext()) {
            //List l=Collections.nCopies(22, rs);
            Symbol symP = (Symbol) symbolsIt.next();
            if (symP.getMatches().toString() == null) { // if it has no matching, it is a unknown symbol "X"
                if (verbose) {
                    System.out.println("[INFO] An unknown symbol was detected in Protein:" + thProtein.getProteinID() + "\t At position:" + symbolIndex);
                }
                if(!symbols.containsKey(symP)) {
                    symbols.put(symP, new LinkedList<Integer>());
                }
                symbols.get(symP).addFirst(symbolIndex);
            }
            symbolIndex++;
        }
        return symbols;
    }

    /**
     * Given a directory of files, check them and generate a report
     *
     * @param dir the directory
     * @return a Report.
     * @throws FileNotFoundException
     * @throws IOException
     */
    public ReportModel getReport(File dir, int missedCleavages) throws FileNotFoundException, IOException {

        int count=0;
        ReportModel rpModel = new ReportModel();
        ProteinProcessor pp = new ProteinProcessor();

        int number_of_proteins = 0;
        long size_of_proteins = 0;
        int number_of_X = 0;
        double maxProteinMass = -1;
        double maxPeptideMass = -1;
        boolean fixedModification=true;
       
        File fList[] = dir.listFiles();
        Random r= new Random();
        
        String folder= "/Users/Eman/Documents/My_Thesis/Thesis-DataSet/ExperimentalMasses/RandomWithCont/";

        for (int i = 0; i < fList.length; i++) {

            if (fList[i].getName().equals(".DS_Store")) {
                continue;
            }
            System.out.println("Loading Sequence file:" + fList[i].getAbsolutePath());
            BufferedReader br = null;

            try {
                br = new BufferedReader(new FileReader(fList[i]));
                Alphabet alpha = AlphabetManager.alphabetForName("PROTEIN");
                SimpleNamespace ns = new SimpleNamespace("biojava");

                RichSequenceIterator iterator = RichSequence.IOTools.readFasta(br, alpha.getTokenization("token"), ns);
                while (iterator.hasNext()) {
                    // in case of more protein sequences in the same file we create a protein for each sequence.
                    // create a new theortical protien
                    RichSequence sequence = iterator.nextRichSequence();
                    Protein protein = new Protein();
                    protein.setMassesList(null);
                    protein.ProteinID = sequence.getAccession();
                    protein.Name= sequence.getName();
                    protein.setFile(fList[i]);
                    protein.setSequence(sequence);
                    protein.setType(Protein.ProteinType.THEORETICAL);
                    //System.out.println("protein Name:"+ protein.getSequence().getName());

                    number_of_proteins++;
                    size_of_proteins += getProteinSize(protein);
                    number_of_X += getUnknowSymbols(protein, false).size();

                    // If the sequence contains U symbol, replace it with A 
                  ArrayList<Integer> uPositions = pp.getPositionsOfU(protein);
                    if (!uPositions.isEmpty()) {
                        protein.setSequence(pp.replacementUsWithOneCode(protein, uPositions));
                    }
                    
                    // If the sequence contains U symbol, replace it with A 
                   ArrayList<Integer> oPositions = pp.getPositionsOfO(protein);
                    if (!oPositions.isEmpty()) {
                        protein.setSequence(pp.replacementOsWithOneCode(protein, oPositions));
                    }
                                      
                    double proteinMass = pp.calcProteinMass(protein);                  
                    pp.digestOneProtein(protein,fixedModification,missedCleavages);
                    double peptideMaxMass = Collections.max(protein.getMassesList());
                    // keep the maximum
                    if (proteinMass > maxProteinMass) {
                        maxProteinMass = proteinMass;
                    }
                    if (peptideMaxMass > maxPeptideMass) {
                        maxPeptideMass = peptideMaxMass;
                    }

                    // to select proteins randomly from the databse and simulate the Mass Spectrometry contaminiation  process////////
                    
//                    if(count>20) System.exit(0);
//                     if(r.nextFloat()<0.0005){
//                                //
//                            double addProb=0.10;
//                            for(float remProb= 0.7f; remProb<=0.91f; remProb+=0.1f){
//                             ArrayList<Double> massListCopy= new ArrayList<Double>(protein.getMassesList().size()+1);
//                             massListCopy.addAll(protein.getMassesList());
//                             
//                                System.out.println("***********************\nRandom Protein:"+ protein.getName() + "\tProtein Mass: "+ protein.getProteinMass() + "\tAdding probability:"+ addProb+ "\tRemving probability:"+ remProb );
//                                writeRandomMassesToFile(massListCopy, folder+protein.Name);
//                                randomize(massListCopy,addProb, remProb); 
//                                writeRandomMassesToFile(massListCopy, folder+protein.Name+ "_"+ protein.getProteinMass() + "_"+ String.valueOf(remProb)+ "_"+ String.valueOf(addProb));
//                              
//                                
//                            }
//                            count++;
//                     }

                    
                    // 1- check how many Xs it has
                    ArrayList<Integer> xPositions = pp.getPositionsOfX(protein);
                    // if is 1 x create a list of proteins and process them instead of the protein
                    if (!xPositions.isEmpty() && xPositions.size() == 1) {
                        ArrayList<Symbol> listofSymbols = pp.getProteinSymbol();

                        for (int j = 0; j < listofSymbols.size(); j++) {
                            // generate a new sequence with X replaced by a symbol form the alphabet 
                            // System.out.println("Symbol:"+ listofSymbols.get(i).getName());
                            Edit edRest = null;
                            try {
                                edRest = new Edit(xPositions.get(0), protein.getSequence().getAlphabet(), listofSymbols.get(j));
                            } catch (IllegalSymbolException ex) {
                                System.err.println(ex);
                            }
                            Protein newProtein = null;
                            try {
                                newProtein = new Protein(protein);
                            } catch (BioException ex) {
                                System.err.println(ex);

                            }
                            try {
                                newProtein.getSequence().edit(edRest);
                            } catch (IndexOutOfBoundsException ex) {
                                System.err.println(ex);
                            } catch (IllegalAlphabetException ex) {
                                System.err.println(ex);
                            } catch (ChangeVetoException ex) {
                                System.err.println(ex);
                            }
                            // calculate all for the new protein
                            proteinMass = pp.calcProteinMass(newProtein);
                            pp.digestOneProtein(newProtein, fixedModification, missedCleavages);
                            peptideMaxMass = Collections.max(newProtein.getMassesList());                           
                            // keep the maximum
                            if (proteinMass > maxProteinMass) {
                                maxProteinMass = proteinMass;
                            }
                            if (peptideMaxMass > maxPeptideMass) {
                                maxPeptideMass = peptideMaxMass;
                            }
                        }
                    }
                }

            } catch (FileNotFoundException ex) {
                // Logger.getLogger(DataReader.class.getName()).log(Level.SEVERE, null, ex);
                System.err.println("File not found!");
            } catch (BioException ex) {
                //Logger.getLogger(DataReader.class.getName()).log(Level.SEVERE, null, ex);
                System.err.println("Can not load sequence Protein Checker:" + ex.getMessage());
            }
        }

        rpModel.setProteinsNO(number_of_proteins);
        rpModel.setXNO(number_of_X);
        rpModel.setdBsize(size_of_proteins);
        rpModel.setMaxProteinMass(maxProteinMass);
        rpModel.setMaxPiptideMass(maxPeptideMass);

        return rpModel;
    }
 
    /**
     * Get the number of unknown symbol "The Xs" in one protein
     *
     * @param thProtein
     * @return Number of Xs in the protein sequence
     */
    public int getNumberOfUnkownSymbolsInProtein(Protein thProtein) {

        int symbolCount = 1;
        RichSequence sequence = thProtein.getSequence();
        Iterator seqSymbols = sequence.iterator();
        while (seqSymbols.hasNext()) {
            Symbol symbol = (Symbol) seqSymbols.next();
            if (symbol.getMatches().toString() == null) {
                symbolCount++;
            }
        }
        return symbolCount;
    }

    /**
     * Get the protein size: the number of characters in the sequence multiply
     * by the 2 "size of char in java 2 bytes"
     *
     * @param thProtein
     * @return
     */
    public int getProteinSize(Protein thProtein) {

        int proteinSize = (thProtein.getSequence().toList().size() * 2);
        return proteinSize;
    }
    

    /**
     * Get the database size: all the characters in the database proteins
     *
     * @param thProteins
     * @return
     */
    public int getDataBaseSize(ArrayList<Protein> thProteins) {
        int dbSize = 0;
        Iterator dbProteinIt = thProteins.iterator();
        while (dbProteinIt.hasNext()) {
            Protein dbP = (Protein) dbProteinIt.next();
            dbSize += getProteinSize(dbP);
        }
        return dbSize;
    }
    
    public int getDBlEngth(ArrayList<Protein> thProteins){
      
        int dbLength=0;
        Iterator dbProteinIt=thProteins.iterator();
        while(dbProteinIt.hasNext()){
            dbLength++;
        }
        
        return dbLength;
    }
    
   
    /**
     * Simulate the noise addition and data loss
     * @param list
     * @param addProbability
     * @param removeProbability
     */
    public void randomize(ArrayList<Double> list, double addProbability, double removeProbability){
        
        // remving 
        int size= list.size();
        LinkedList<Double> rem= new LinkedList<Double> ();
        Random r= new Random();
        for(int i=0;i<size;i++){
            
        
            if(r.nextFloat()<removeProbability) {
                rem.addFirst(list.get(i));
            }
        }
        
        list.removeAll(rem);
        
        
        Set<Double> add= new LinkedHashSet<Double> ();
        for(int i=0;i<size;i++){
         
            if(r.nextFloat()<addProbability){
                double randomToAdd=(400 + (7000 - 400) * r.nextFloat());
                    
                        add.add(randomToAdd);    
            }
        }
        
      list.addAll(add);
    }
    
    public void writeRandomMassesToFile(ArrayList<Double> randomMasses, String randomMassesfile) {
        PrintWriter pw = null;
        try {
            pw = new PrintWriter(new FileOutputStream(randomMassesfile));
            Iterator it = randomMasses.iterator();
            while (it.hasNext()) {
                double randomMass = (Double) it.next();
                pw.println(randomMass);
            }
            pw.flush();
            pw.close();
        } catch (FileNotFoundException ex) {
            Logger.getLogger(ProteinChecker.class.getName()).log(Level.SEVERE, null, ex);
        } finally {
            pw.close();
        }
    }
    
    
}
