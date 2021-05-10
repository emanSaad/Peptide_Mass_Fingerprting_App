/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package org.ualg.pt.pmfgui;

import java.awt.Color;
import java.io.BufferedReader;
import java.io.File;
import java.io.FileNotFoundException;
import java.io.FileReader;
import java.io.IOException;
import java.text.SimpleDateFormat;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Calendar;
import java.util.Iterator;
import java.util.LinkedHashSet;
import java.util.LinkedList;
import java.util.Map;
import java.util.Set;
import java.util.logging.Level;
import java.util.logging.Logger;
import javax.swing.JOptionPane;
import org.biojava.bio.BioException;
import org.biojava.bio.symbol.Alphabet;
import org.biojava.bio.symbol.AlphabetManager;
import org.biojava.bio.symbol.Symbol;
import org.biojavax.SimpleNamespace;
import org.biojavax.bio.seq.RichSequence;
import org.biojavax.bio.seq.RichSequenceIterator;
import org.netbeans.api.progress.ProgressHandle;
import org.netbeans.api.progress.ProgressHandleFactory;
import org.openide.util.Cancellable;
import org.openide.util.Exceptions;
import org.openide.util.RequestProcessor;
import org.openide.util.TaskListener;
import org.ualg.pt.pmfambiguity.AmbiguityResolver;
import org.ualg.pt.pmfcore.DataReader;
import org.ualg.pt.pmfcore.Match;
import org.ualg.pt.pmfcore.MatchComparator;
import org.ualg.pt.pmfcore.Protein;
import org.ualg.pt.pmfcore.ProteinChecker;
import org.ualg.pt.pmfcore.ProteinMatcher;
import org.ualg.pt.pmfcore.ProteinProcessor;
import org.ualg.pt.pmfcore.ReportModel;

/**
 *
 * @author Eman
 */
public abstract class ProteinFinder {

    private final static RequestProcessor RP = new RequestProcessor("interruptible tasks", 1, true);
    static RequestProcessor.Task theTask = null;
    static Runnable runnable = null;

    /**
     * @param args the command line arguments
     */
    public static void search( final MainTopComponent parent) throws FileNotFoundException, IOException, BioException {

        // progressbar handelr in the GUI
        final ProgressHandle ph = ProgressHandleFactory.createHandle("Collecting information ...", new Cancellable() {

            @Override
            public boolean cancel() {

                if (null == theTask) {
                    return false;
                }

                return theTask.cancel();

            }
        });

        // TODO code application logic here
        final ProteinProcessor proteinprocessor = new ProteinProcessor();
        final File proteinDB = parent.getProteinDB();
        final File decoyDB = parent.getDecoyDB();
        final File massesDir = parent.getMassList();
        final File contaminantDB = parent.getContaminantsList();
        //File massesFile = new File("/home/faroq/Documents/Eman/SwissSamplesForTest/Q6GZV2");

        final File[] massesFiles = massesDir.listFiles();

        String timeStamp = new SimpleDateFormat("yyyy/MM/dd HH:mm:ss").format(Calendar.getInstance().getTime());
        System.out.println("Time of Today: " + timeStamp);
        parent.appendResult("Date:" + timeStamp + "\n", Color.black);

        // this has to be option in the software interface
        final double contaminantsThreshold = parent.getContaminantsThreshold();
        final double dataBasethreshold = parent.getDatabaseThreshold(); // matching threshold in  Da
        final int missedCleavages = parent.getNoMissedCleavages();
        final boolean fixedModification = parent.isFixedModified();
        final boolean variableModification = parent.isVariableModified();
        final int reportHits = parent.getMaximumNoHits();
        final float significanceThreshold = parent.getSignificanceThreshold();

        System.out.println("Settings:");
        System.out.println("Contaminant Threshold: " + contaminantsThreshold);
        System.out.println("DataBase Threshold: " + dataBasethreshold + "Da");
        System.out.println("Number of Missed Cleavages: " + missedCleavages);
        System.out.println("Fixed Modification: " + fixedModification);
        System.out.println("Variable Modification: " + variableModification);
        System.out.println("Number of Top Hits:" + reportHits);
        System.out.println("Significance threshold:" + significanceThreshold);

        parent.appendResult("Settings:\n", Color.black);
        parent.appendResult("\t-Protein DB:" + proteinDB.getAbsolutePath() + "\n", Color.gray);
        parent.appendResult("\t-Decoy DB:" + decoyDB.getAbsolutePath() + "\n", Color.gray);
        parent.appendResult("\t-Mass List Folder:" + massesDir.getAbsolutePath() + "\n", Color.gray);
        parent.appendResult("\t-Contaminant Threshold: " + contaminantsThreshold + "\n", Color.gray);
        parent.appendResult("\t-DataBase Threshold: " + dataBasethreshold + "Da\n", Color.gray);
        parent.appendResult("\t-Number of Missed Cleavages: " + missedCleavages + "\n", Color.gray);
        parent.appendResult("\t-Fixed Modification: " + fixedModification + "\n", Color.gray);
        parent.appendResult("\t-Variable Modification: " + variableModification + "\n", Color.gray);
        parent.appendResult("\t-Number of Top Hits:" + reportHits + "\n", Color.gray);
        parent.appendResult("\t-Significance threshold:" + significanceThreshold + "\n", Color.gray);
        
        ph.start();
        ph.switchToDeterminate(proteinDB.listFiles().length + decoyDB.listFiles().length + massesDir.listFiles().length);
      

        // put them in a thread
        runnable = new Runnable() {
            @Override
            public void run() {
                int workUnitSofar = proteinDB.listFiles().length;
                // check the sequence database and collect some needed information
                System.out.println("Checking protien database ...");
                parent.appendResult("Checking protien database ... \n", Color.black);
                ph.setDisplayName("Checking protien database ... ");
                ProteinChecker proteinChecker = new ProteinChecker();
                ReportModel rpmModel = null;
                try {
                    rpmModel = proteinChecker.getReport(proteinDB, missedCleavages);
                } catch (IOException ex) {
                    Exceptions.printStackTrace(ex);
                }
                System.out.println("In \tUnknown symbols: " + rpmModel.getXNO());

                ph.progress(workUnitSofar);

                System.out.println("Checking decoy database ...");
                parent.appendResult("Checking decoy database ... \n", Color.black);
                ph.setDisplayName("Checking decoy database ...");
                ProteinChecker proteinChecker1 = new ProteinChecker();
                ReportModel rpmModel1 = null;
                try {
                    rpmModel1 = proteinChecker1.getReport(decoyDB, missedCleavages);
                } catch (IOException ex) {
                    Exceptions.printStackTrace(ex);
                }
                System.out.println("In \tUnknown symbols: " + rpmModel1.getXNO());
                workUnitSofar = workUnitSofar + decoyDB.listFiles().length;
                ph.progress(workUnitSofar);

                for (File massesFile : massesFiles) {
                    try {
                        // read the masses
                        if (Thread.interrupted()) {
                            JOptionPane.showMessageDialog(null, "Search process was interrupted");
                            return;
                        }

                        System.out.println("Loaded masses...... from file:" + massesFile.getPath());
                        ph.setDisplayName("Searching for " + massesFile.getName() + " ...");
                        parent.appendResult("Searching for " + massesFile.getName() + " ...", Color.BLUE, Color.lightGray);
                        parent.appendResult(" ", Color.black, Color.WHITE);
                        Protein massFile = null;
                        try {
                            massFile = proteinprocessor.readMassesFile(massesFile);
                        } catch (IOException ex) {
                            Exceptions.printStackTrace(ex);
                        }
                        ArrayList<Protein> massesList = new ArrayList<Protein>();
                        massesList.add(massFile);

                        // collect significace
                        System.out.println("Matching with database");
                        Match[] swissMatches = ProteinFinder.batchProcess(proteinDB, rpmModel, massesList, contaminantDB, contaminantsThreshold, dataBasethreshold, fixedModification, variableModification, missedCleavages);
                        System.out.println("Matching with decoy");
                        Match[] randomSwissMatches = ProteinFinder.batchProcess(decoyDB, rpmModel, massesList, contaminantDB, contaminantsThreshold, dataBasethreshold, fixedModification, variableModification, missedCleavages);

                        Arrays.sort(swissMatches, new MatchComparator());
                        Arrays.sort(randomSwissMatches, new MatchComparator());
                        //SignificanceCalculator significanceCalculator = new SignificanceCalculator();
                        SignificanceCalculator.calcSignificance(swissMatches, randomSwissMatches, reportHits, 1 - significanceThreshold, parent);
                        workUnitSofar = workUnitSofar + 1;
                        ph.progress(workUnitSofar);
                        System.gc(); //  force the cleaning of any used memory by previouse search

                    } catch (IOException ex) {
                        Exceptions.printStackTrace(ex);
                    }

                }
            }
        };

        theTask = RP.create(runnable);
        theTask.addTaskListener(new TaskListener() {
            @Override
            public void taskFinished(org.openide.util.Task task) {
                //make sure that we get rid of the ProgressHandle
                //when the task is finished
                ph.finish();
                JOptionPane.showMessageDialog(null, "Searching has finished!");

            }
        });

        //start the progresshandle the progress UI will show 500s after
          //this actually start the task
        theTask.schedule(0);

    }

    public static Match[] batchProcess(File dbdir, ReportModel rpModel, ArrayList<Protein> massesList, File contaminantsdir, double contaminantsThreshold, double dataBasethreshold, boolean fixedModification, boolean variableModification, int missedCleavages) throws IOException {
        ProteinProcessor pp = new ProteinProcessor();

        /**
         * *********************************** first remove all contaminants
         * ************************************
         */
        System.out.println("Removing contaminants ...................");
        ArrayList<Double> correctMasslist = new ArrayList<Double>();
        File fList[] = contaminantsdir.listFiles();
        int dblength = 0;
        // compare with all masses
        Iterator exIterator = massesList.iterator();
        while (exIterator.hasNext()) {
            Protein exP = (Protein) exIterator.next(); // get the current protein
            Iterator massExIterator = exP.getMassesList().iterator(); // get masses list of the current protein
            while (massExIterator.hasNext()) {
                // Get current mass then compare it to all other masses in the other proteins
                double currentExMass = (Double) massExIterator.next();
                for (int i = 0; i < fList.length; i++) {

                    if (fList[i].getName().equals(".DS_Store")) {
                        continue;
                    }
                    // System.out.println("Loading Sequence file:" + fList[i].getAbsolutePath());

                    RichSequence rec = null;
                    BufferedReader br = null;
                    try {
                        br = new BufferedReader(new FileReader(fList[i]));
                        Alphabet alpha = AlphabetManager.alphabetForName("PROTEIN");
                        SimpleNamespace ns = new SimpleNamespace("biojava");
                        Protein newExProtein = new Protein();
                        RichSequenceIterator iterator = RichSequence.IOTools.readFasta(br, alpha.getTokenization("token"), ns);

                        while (iterator.hasNext()) {
                            // create a new theortical protien
                            ProteinProcessor proteinprocessor = new ProteinProcessor();
                            RichSequence sequence = iterator.nextRichSequence();
                            Protein protein = new Protein();
                            protein.setMassesList(null);
                            protein.setProteinID(sequence.getName());
                            protein.setName(sequence.getAccession());
                            protein.setFile(fList[i]);
                            protein.setSequence(sequence);
                            protein.setType(Protein.ProteinType.EXPERIMENTAL);

                            pp.digestOneProtein(protein, false, 0);
                            //print contaminants masses list.
                            //System.out.println("Contaminants proteins masses list: " + protein.getMassesList());
                            Iterator massContIterator = protein.getMassesList().iterator();
                            while (massContIterator.hasNext()) {
                                double currentcontMass = (Double) massContIterator.next();
                                // Do the comparision
                                double result = Math.abs(currentcontMass - currentExMass);
                                if (result <= contaminantsThreshold) {
                                    // add the experimental masses that are appear to be contaminants to contaminants list
                                    correctMasslist.add(currentExMass);
                                    //exP.getMassesList().remove(currentExMass);                                   
                                }
                            }
                        }

                    } catch (FileNotFoundException ex) {
                        Logger.getLogger(DataReader.class.getName()).log(Level.SEVERE, null, ex);
                        System.err.println("File not found!");
                    } catch (BioException ex) {
                        Logger.getLogger(DataReader.class.getName()).log(Level.SEVERE, null, ex);
                        System.err.println("Can not load contaminants sequence:" + ex.getMessage());
                    }
                }
            }

            //remove the masses that are contaminants, and add the correct experimental masses to the massesList to match them with theoretical proteins
            // System.out.println("experimental masses before removing: " + exP.getMassesList());
            exP.getMassesList().removeAll(correctMasslist);
            System.out.println("valid masses List: " + exP.getMassesList());

        }

        /**
         * *********************************** build the frequency matrix
         * ***********************************
         */
        System.out.println("Building the frequencey matrix .....");
        ProteinMatcher proteinMatcher = new ProteinMatcher((int) Math.round(rpModel.getMaxProteinMass()), (int) Math.round(rpModel.getMaxPiptideMass()));

        fList = dbdir.listFiles();

        for (int i = 0; i < fList.length; i++) {

            if (fList[i].getName().equals(".DS_Store")) {
                continue;
            }
            System.out.println("Loading Sequence file:" + fList[i].getAbsolutePath());
            RichSequence rec = null;
            BufferedReader br = null;
            try {
                br = new BufferedReader(new FileReader(fList[i]));
                Alphabet alpha = AlphabetManager.alphabetForName("PROTEIN");
                SimpleNamespace ns = new SimpleNamespace("biojava");
                RichSequenceIterator iterator = RichSequence.IOTools.readFasta(br, alpha.getTokenization("token"), ns);
                while (iterator.hasNext()) {
                    // create a new theortical protien
                    RichSequence sequence = iterator.nextRichSequence();
                    Protein protein = new Protein();
                    protein.setMassesList(null);
                    protein.setProteinID(sequence.getAccession());
                    protein.setFile(fList[i]);
                    protein.setSequence(sequence);
                    protein.setType(Protein.ProteinType.THEORETICAL);
                    protein.setName(sequence.getAccession());

                    //count the number of proteins in database
                    dblength++;
                    // replace any unknow symbols
                    AmbiguityResolver ambiguityResolver = new AmbiguityResolver();
                    Map<Symbol, LinkedList<Integer>> unknowns = ProteinChecker.getUnknowSymbolsMap(protein, false);

                    boolean hasSingleUnknowns = false;
                    int counter = 0;
                    Set<Double> peptideMasses = new LinkedHashSet<Double>();
                    double averageMass = 0.0;
                    double numberOfReplacements = 0;
                    if (!unknowns.isEmpty()) {
                        ambiguityResolver.calcMaps(protein, unknowns);
                        Protein newProtein = ambiguityResolver.resolveSequence(protein);

                        while (newProtein != null) {
                            // System.out.println("Processing Generated Protein: "+ newProtein );
                            // System.out.println("counter:"+ ++counter);
                            numberOfReplacements++;
                            hasSingleUnknowns = true;
                            // calcuate the new masslist 
                            pp.digestOneProtein(newProtein, fixedModification, missedCleavages);
                            pp.calcProteinMass(newProtein);
                            // System.out.println("Protein Digest:"+newProtein.getMassesList());
                            peptideMasses.addAll(newProtein.getMassesList());
                            averageMass += newProtein.getProteinMass();
                            // here we have new protein digested and ready to go the matching matrix

                            // System.out.println("New Protein ID: " +protein.getProteinID()+"\tProtein Mass:"+ newProtein.getProteinMass());
                            //proteinMatcher.printFrequenceyTable();
                            newProtein = ambiguityResolver.resolveSequence(protein);

                        }
                    }
                    if (!hasSingleUnknowns) {

                        // no Xes or more than one X do the natural digestion
                        // 1 - digest
                        pp.digestOneProtein(protein, fixedModification, missedCleavages);
                        pp.calcProteinMass(protein);
                        proteinMatcher.addFrequencies(protein);
                        proteinMatcher.addFrequenciesforProteins(protein);

                        //System.out.println("Protein ID: " +protein.getProteinID()+"\tProtein Mass:"+ protein.getProteinMass());
                    } else { // it has some unknown being replaced,
                        // add the frequencies here
                        // use the average mass of the protein with the collected unique peptide masses.
                        proteinMatcher.addFrequencies(peptideMasses, averageMass / numberOfReplacements);
                        proteinMatcher.addFrequenciesforProteins(protein);
                    }

                }
            } catch (FileNotFoundException ex) {
                Logger.getLogger(DataReader.class.getName()).log(Level.SEVERE, null, ex);
                System.err.println("File not found!");
            } catch (BioException ex) {
                Logger.getLogger(DataReader.class.getName()).log(Level.SEVERE, null, ex);
                System.err.println("Can not load  DB (in add frequency)sequence:" + ex.getMessage());
            }
        }

        proteinMatcher.normalizeToLargestValue();

        /**
         * ***********************************finally do the matching
         * ***************************************
         */
        System.out.println("Matching..................");

        // the kbest protein list
        LinkedList<Match> kbest = new LinkedList<Match>();
        double msMassesLength = 0;

        fList = dbdir.listFiles();
        System.out.println("Loading Sequence file:" + dbdir.getAbsolutePath());
        // compare with all masses
        Iterator exIterator1 = massesList.iterator();
        while (exIterator1.hasNext()) {
            Protein exP = (Protein) exIterator1.next(); // get the current protein
            msMassesLength = exP.getMassesList().size();
            Iterator massExIterator = exP.getMassesList().iterator(); // get masses list of the current protein

            for (int i = 0; i < fList.length; i++) {

                if (fList[i].getName().equals(".DS_Store")) {
                    continue;
                }
                // System.out.println("Loading Sequence file:" + fList[i].getAbsolutePath());

                int j = 0;
                BufferedReader br = null;
                try {
                    br = new BufferedReader(new FileReader(fList[i]));
                    Alphabet alpha = AlphabetManager.alphabetForName("PROTEIN");
                    SimpleNamespace ns = new SimpleNamespace("biojava");

                    RichSequenceIterator iterator = RichSequence.IOTools.readFasta(br, alpha.getTokenization("token"), ns);
                    while (iterator.hasNext()) {
                        // create a new theortical protien
                        RichSequence sequence = iterator.nextRichSequence();
                        Protein protein = new Protein();
                        protein.setMassesList(null);
                        protein.setProteinID(sequence.getName());
                        protein.setFile(fList[i]);
                        protein.setSequence(sequence);
                        protein.setType(Protein.ProteinType.THEORETICAL);
                        protein.setName(sequence.getAccession());

                        // 1- check how many Xs it has
                        AmbiguityResolver ambiguityResolver = new AmbiguityResolver();
                        Map<Symbol, LinkedList<Integer>> unknowns = ProteinChecker.getUnknowSymbolsMap(protein, false);
                        boolean hasSingleUnknowns = false;
                        int counter = 0;
                        if (!unknowns.isEmpty()) {
                            ambiguityResolver.calcMaps(protein, unknowns);
                            Protein newProtein = ambiguityResolver.resolveSequence(protein);

                            while (newProtein != null) {
                                hasSingleUnknowns = true;

                                if (variableModification == true) {

                                    // calcuate the new masslist 
                                    pp.digestionAndCalculateMassesWithModification(newProtein, fixedModification, variableModification, missedCleavages);
                                } else {
                                    pp.digestOneProtein(newProtein, fixedModification, missedCleavages);

                                }
                                // System.out.println("Protein Digest:"+newProtein);
                                pp.calcProteinMass(newProtein);
                                // match with newprotein
                                double score = proteinMatcher.doMatchingOneByOne(exP, newProtein, dataBasethreshold);
                                if (newProtein.getMatchedCoverage() >= 0.18) {
                                    kbest.addLast(new Match(newProtein.getProteinID() + ambiguityResolver.getMapAsString(), score, newProtein, msMassesLength));
                                }

                                //this.insert(newProtein.getProteinID() + ambiguityResolver.getMapAsString(), newProtein, score, msMassesLength, kbest);
                                //}
                                newProtein = ambiguityResolver.resolveSequence(protein);
                                // here we have new protein digested and ready to go the matching matrix   
                            }
                        }

                        if (!hasSingleUnknowns) {

                            // no Xes or more than one X do the natural digestion
                            // 1 - digest
                            if (variableModification == true) {

                                pp.digestionAndCalculateMassesWithModification(protein, fixedModification, variableModification, missedCleavages);
                            } else {
                                pp.digestOneProtein(protein, fixedModification, missedCleavages);
                            }
                            pp.calcProteinMass(protein);

                            // match with the protein
                            double score = proteinMatcher.doMatchingOneByOne(exP, protein, dataBasethreshold);
                            double mass = protein.getProteinMass();
                            if (protein.getMatchedCoverage() >= 0.18) {
                                kbest.addLast(new Match(protein.getProteinID(), score, protein, msMassesLength));
                            }
                            // this.insert(protein.getProteinID(), protein, score, msMassesLength, kbest);
                            //System.out.println(protein.getName()+"\t\t"+protein.getProteinID()+score +" \t\t ");
//                            (protein.getMatchedCoverage()==0.18)
//                                  System.out.println("Protein Score with 0.18 threshold: " + protein.getProteinID()+"\t"+score );
                        }
                    }
                } catch (FileNotFoundException ex) {
                    // Logger.getLogger(DataReader.class.getName()).log(Level.SEVERE, null, ex);
                    System.err.println("File not found!");
                } catch (BioException ex) {
                    //Logger.getLogger(DataReader.class.getName()).log(Level.SEVERE, null, ex);
                    System.err.println("Can not load  DB sequence in matching:" + ex.getMessage());
                }
            }

        }
        int length = kbest.size();
        Match[] matches = new Match[length];
        for (int i = 0; i < length; i++) {
            matches[i] = kbest.get(i);
        }

        return matches;

    }

    public void insert(String proteinID, Protein protein, double proteinScore, double experementals, Match[] kbest) {

        // to be used in the interface for filtering
        // coverage threshold
        double coverageThreshold = 0.18;
        //Edit
        if (protein.getMatchedCoverage() < coverageThreshold) {
            proteinScore = 0;
        }

        //End edit 
        for (int i = 0; i < kbest.length; i++) {

            // two cases
            if (kbest[i] == null) { // 1- empty
                kbest[i] = new Match(proteinID, proteinScore, protein, experementals);
                break;
            }

            // 2- nonempty
            if (proteinScore >= kbest[i].getScore()) {
                if (i == kbest.length - 1) { // if is out of bounds, replace it 
                    kbest[i] = new Match(proteinID, proteinScore, protein, experementals);
                } else {

                    for (int j = kbest.length - 1; j > i; j--) {   // shift all less elements
                        if (kbest[j - 1] == null) {
                            continue;
                        }
                        kbest[j] = kbest[j - 1].copy();
                    }
                    kbest[i] = new Match(proteinID, proteinScore, protein, experementals);
                    break;
                }
            }
        }
    }

}
