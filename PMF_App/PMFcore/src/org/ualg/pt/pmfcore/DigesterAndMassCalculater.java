/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package org.ualg.pt.pmfcore;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashSet;
import java.util.Iterator;
import java.util.LinkedList;
import java.util.List;
import java.util.Set;
import org.apache.commons.lang.ArrayUtils;
import org.biojava.bio.BioException;
import org.biojava.bio.proteomics.Digest;
import org.biojava.bio.proteomics.MassCalc;
import org.biojava.bio.proteomics.ProteaseManager;
import org.biojava.bio.seq.ProteinTools;
import org.biojava.bio.symbol.FiniteAlphabet;
import org.biojava.bio.symbol.IllegalSymbolException;
import org.biojava.bio.symbol.Symbol;
import org.biojava.bio.symbol.SymbolList;
import org.biojava.bio.symbol.SymbolPropertyTable;
import org.biojava.utils.ChangeVetoException;
import org.biojavax.bio.seq.RichSequence;
import org.biojavax.bio.seq.SimpleRichFeature;

/**
 *
 * @author Eman
 */
public class DigesterAndMassCalculater {

    /**
     * Pass an a protein sequence, calculate the mass value for the protein
     * (molecule weight)
     *
     * @param sequence
     * @return molecule weight of protein sequence
     * @throws IllegalSymbolException
     */
    public static double getProteinMass(RichSequence sequence) throws IllegalSymbolException {

        MassCalculator massC = new MassCalculator(SymbolPropertyTable.MONO_MASS, true);  
        double pMass = massC.getMolecularWeight2(sequence);
        return pMass;
    }
    
    /**
     *Get the symbols of ambiguous symbol B which are D and N and put them in arrayList
     * @return Corresponding symbols of B
     */
    public static Set<Symbol> getBsymbols() {
        Set <Symbol> symbolsOfB = new HashSet<Symbol>();
        FiniteAlphabet protein = ProteinTools.getAlphabet();
        Iterator proteinSym = protein.iterator();
        while (proteinSym.hasNext()) {
            Symbol s = (Symbol) proteinSym.next();
            if (s.getName().equals("ASN") || s.getName().equals("ASP")) {
                symbolsOfB.add(s);
                //System.out.println("List of protein symbols: " + s.getName());
            }
        }
        return symbolsOfB;
    }
    
    public static Set<Symbol> getXsymbols(){
        Set <Symbol> symbolsOfX = new HashSet<Symbol>();
        FiniteAlphabet protein = ProteinTools.getAlphabet();
        Iterator proteinSym = protein.iterator();
        while (proteinSym.hasNext()) {
            Symbol s = (Symbol) proteinSym.next();
                symbolsOfX.add(s);
                //System.out.println("List of protein symbols: " + s.getName());      
        }
        return symbolsOfX;       
    }
    
    /**
     *Get the symbols of ambiguous symbol Z which are E and Q and put them in arrayList
     * @return Corresponding symbols of Z
     */
    public static Set<Symbol> getZsymbols() {
        Set <Symbol> symbolsOfZ = new HashSet<Symbol>();
        FiniteAlphabet protein = ProteinTools.getAlphabet();
        Iterator proteinSym = protein.iterator();
        while (proteinSym.hasNext()) {
            Symbol s = (Symbol) proteinSym.next();
            if (s.getName().equals("GLU") || s.getName().equals("GLN")) {
                symbolsOfZ.add(s);
                //System.out.println("List of protein symbols: " + s.getName());
            }
        }
        return symbolsOfZ;
    }
    
    /**
     *Get the symbols of ambiguous symbol J which are I and L and put them in arrayList
     * @return Corresponding Symbols of J
     */
    public static Set<Symbol> getJsymbols() {
        Set <Symbol> symbolsOfJ = new HashSet<Symbol>();
        FiniteAlphabet protein = ProteinTools.getAlphabet();
        Iterator proteinSym = protein.iterator();
        while (proteinSym.hasNext()) {
            Symbol s = (Symbol) proteinSym.next();
            if (s.getName().equals("ILE") || s.getName().equals("LEU")) {
                symbolsOfJ.add(s);
                //System.out.println("List of protein symbols: " + s.getName());
            }
        }
        return symbolsOfJ;
    }
    
    /**
     * Digest the proteins to peptides and calculate their masses.
     *
     * @param seq
     * @param enzymeType
     * @return
     * @throws BioException
     * @throws ChangeVetoException
     */
    public static ArrayList<Double> proteinDigestionAndMassesCalculation1(Protein sourceprotein, String enzymeType,boolean fixedModification, int maxMissedCleavages) throws BioException, ChangeVetoException {
        // here 0 is the number of missed cleaveges sites, when build the interface make setMaxMissedCleavages to the user option.
        //int maxMissedCleavages = 1;
        RichSequence seq=sourceprotein.getSequence();
        Digest digest = new Digest();
        digest.setMaxMissedCleavages(maxMissedCleavages);
        digest.setProtease(ProteaseManager.getProteaseByName(enzymeType));
        digest.setSequence(seq);
        digest.addDigestFeatures();
        ArrayList<Double> thMasses = new ArrayList<Double>();
        Iterator iterator = digest.getSequence().features();
        
        LinkedList<Peptide> peptides= new LinkedList<Peptide>();
        
        FiniteAlphabet protein = ProteinTools.getAlphabet();       
        double mass = 0d;
         
         MassCalc mc = new MassCalc(SymbolPropertyTable.MONO_MASS, true);
        while (iterator.hasNext()) {
            // iterate through the features
            SimpleRichFeature f = (SimpleRichFeature) iterator.next();
            SymbolList symList = f.getSymbols();
            //get pI value for the peptide
          
            Iterator peptideIt = symList.iterator();
            while (peptideIt.hasNext()) {
                Symbol oneCode = (Symbol) peptideIt.next();
                            
                mc.setSymbolModification('X', 110);
                mc.setSymbolModification('B', 132.6108);                 
                mc.setSymbolModification('Z', 146.6375);       
                mc.setSymbolModification('J', 131.1736);
                mc.setSymbolModification('U', 131.1736);
                mc.setSymbolModification('O', 100.1736);
                
                if (oneCode.getName().equals("CYS") && fixedModification==true)
                {         
                    //System.out.println("Here is fixed modification......" + oneCode.getName().toString());
                    //mc.setSymbolModification('C', 175.159);
                      mc.setSymbolModification('C', 160.04323);
                }                
            }
                     
            mass = mc.getMass(symList);           
            thMasses.add(mass);
           //System.out.println("Position: "+f.getLocation() +"\t"+f.getSymbols().seqString()+ "\t"+"Symbol at position two: "+f.getLocation().contains(3) +"\t"+mass);  
            // System.out.println("Position: " + f.getLocation() + " " + symList.seqString() +"   " + mass + "\n");

            Peptide peptide=new  Peptide();
                 peptide.peptideValue=f.getSymbols().seqString();
                 peptide.mass=mass;
                 peptide.start=f.getLocation().getMin();
                 peptide.end=f.getLocation().getMax();
                 
                 peptides.addLast(peptide);
        }
        sourceprotein.setPeptides(peptides);
        //System.out.println("Return All masses: " + thMasses+"Protein Name: " + seq.getName().toString());

       // System.out.println("All Masses: " + thMasses);
        return thMasses;
        
    }

    
    /**
     * Get the peptides of digested sequence
     * @param seq
     * @param enzymeType
     * @param fixedModification
     * @param maxMissedCleavages
     * @return
     * @throws BioException
     * @throws ChangeVetoException
     */
    public static LinkedList<Peptide> getPeptides( RichSequence seq, String enzymeType,boolean fixedModification,boolean variableModification, int maxMissedCleavages) throws BioException, ChangeVetoException {
        
        
        Digest digest = new Digest();
        digest.setMaxMissedCleavages(maxMissedCleavages);
        digest.setProtease(ProteaseManager.getProteaseByName(enzymeType));
        digest.setSequence(seq);
        digest.addDigestFeatures();
        ArrayList<Double> thMasses = new ArrayList<Double>();
        Iterator iterator = digest.getSequence().features();

        LinkedList<Peptide> peptides = new LinkedList<Peptide>();
        FiniteAlphabet protein = ProteinTools.getAlphabet();
        double pI = 0d;
        double mass = 0d;
        double[] vMasses= new double[1];
        double[] newVariableMasses;
        MassCalc mc = new MassCalc(SymbolPropertyTable.MONO_MASS, true);
        while (iterator.hasNext()) {
            // iterate through the features
            SimpleRichFeature f = (SimpleRichFeature) iterator.next();
            SymbolList symList = f.getSymbols();

            Iterator peptideIt = symList.iterator();
            while (peptideIt.hasNext()) {
                Symbol oneCode = (Symbol) peptideIt.next();
                mc.setSymbolModification('X', 110);
                mc.setSymbolModification('B', 132.6108);

                mc.setSymbolModification('Z', 146.6375);

                mc.setSymbolModification('J', 131.1736);

                mc.setSymbolModification('U', 131.1736);

                mc.setSymbolModification('O', 100.1736);

                if (oneCode.getName().equals("CYS") &&fixedModification == true) {
                   // System.out.println("Fixed modification in peptides... " + oneCode.getName().toString());
                    mc.setSymbolModification('C', 160.04323);
                }
                if(oneCode.getName().equals("Met") && variableModification==true){
                    mc.addVariableModification('M', vMasses);
                }
            }

            //mass = mc.getMass(symList);
            newVariableMasses = mc.getVariableMasses(symList);            
             for(int i=0; i<newVariableMasses.length;i++) {
                 thMasses.add(newVariableMasses[i]);
          //  System.out.println("Position: " + f.getLocation() + " " + symList.seqString() +"   " + i +"  "+ newVariableMasses[i] + "\n");
            }
            //thMasses.add(mass);
            //System.out.println("Position: "+f.getLocation() +"\t"+f.getSymbols().seqString()+ "\t"+"Symbol at position two: "+f.getLocation().contains(3) +"\t"+mass);  
            Peptide peptide = new Peptide();
            peptide.peptideValue = f.getSymbols().seqString();
            peptide.mass = mass;
            peptide.start = f.getLocation().getMin();
            peptide.end = f.getLocation().getMax();
            peptides.addLast(peptide);
        }
        return peptides;
    }

                 
    /**
     * Digest the proteins to peptides and calculate their masses when modifications are selected.
     *
     * @param seq
     * @param enzymeType
     * @return
     * @throws BioException
     * @throws ChangeVetoException
     */
    public static ArrayList<Double> proteinDigestionAndMassesCalculationWithModification(Protein sourceprotein, String enzymeType, boolean fixedModificationSelected, boolean variableModificationSelected, int maxMissedCleavages) throws BioException, ChangeVetoException {
        // here 0 is the number of missed cleaveges sites, when build the interface make setMaxMissedCleavages to the user option.
       RichSequence seq=sourceprotein.getSequence();
        //double vOldMass=149.2124;
        //double vNewMass=149.2124+16;
        //double[] variableMasses = {149.2124, 165.2124};
        double [] variableMasses=new double[1];
        variableMasses[0]=131.04049+16;
        Symbol symbolB;
        Symbol symbolZ;
        Symbol symbolJ;
        double mass=0;
        
        List<Double> massesArrayList=new ArrayList<Double>();
        ArrayList<Double> thMasses = new ArrayList<Double>();
        Set <Symbol> symbolsOfB=getBsymbols();
        Set <Symbol> symbolsOfZ=getZsymbols();
        Set <Symbol> symbolsOfJ=getJsymbols();
        double[] newVariableMasses;
       // List <Double> massesWhenModification=Arrays.asList(newVariableMasses);
        
        FiniteAlphabet protein = ProteinTools.getAlphabet();
        Digest digest = new Digest();
        digest.setMaxMissedCleavages(maxMissedCleavages);
        digest.setProtease(ProteaseManager.getProteaseByName(enzymeType));
        digest.setSequence(seq);
        digest.addDigestFeatures();
         LinkedList<Peptide> peptides= new LinkedList<Peptide>();
        Iterator iterator = digest.getSequence().features();
        while (iterator.hasNext()) {
            // iterate through the features
            SimpleRichFeature f = (SimpleRichFeature) iterator.next();
            SymbolList symList = f.getSymbols();
            MassCalc mc = new MassCalc(SymbolPropertyTable.MONO_MASS, true);
            Iterator peptideIt = symList.iterator();
            while (peptideIt.hasNext()) {
                Symbol oneCode = (Symbol) peptideIt.next();         
                    mc.setSymbolModification('X', 110);             
                    mc.setSymbolModification('B', 132.6108);                                 
                    mc.setSymbolModification('Z', 146.6375);                 
                    mc.setSymbolModification('J', 131.1736);
                    mc.setSymbolModification('U', 131.1736);
                    mc.setSymbolModification('O', 100.1736);
    
                if (oneCode.getName().equals("CYS") && fixedModificationSelected==true) {                   
                    mc.setSymbolModification('C', 160.04323);
                    //System.out.println("Fixed modification in digestion with modification method.");
                }                          
            }
            
            //mass = mc.getMass(symList);
            //System.out.println("Symbol list: "+symList.seqString()+"  masses: " + mass);
           // thMasses.add(mass);
            
            if(variableModificationSelected==true)
            {
                
                //variableMasses[1]=vNewMass;
               // mc.setSymbolModification('M', variableMasses[0]);
                //for(int i=0;i<variableMasses.length;i++)
                //{
                  //variableMasses[0]=vOldMass;
                  
                //}
                
                mc.addVariableModification('M', variableMasses);
            }
            
            newVariableMasses = mc.getVariableMasses(symList);            
             for(int i=0; i<newVariableMasses.length;i++) {
                 thMasses.add(newVariableMasses[i]);
          //  System.out.println("Position: " + f.getLocation() + " " + symList.seqString() +"   " + i +"  "+ newVariableMasses[i] + "\n");
            }
            //***********************
        Peptide peptide=new  Peptide();
                 peptide.peptideValue=f.getSymbols().seqString();
                 peptide.mass=mass;
                 peptide.start=f.getLocation().getMin();
                 peptide.end=f.getLocation().getMax();            
                 peptides.addLast(peptide);
        //***********************
        } 
         sourceprotein.setPeptides(peptides);
       // System.out.println("Return All masses: " + thMasses+"Protein Name: " + seq.getName().toString());
      
        return thMasses;
    }   
}
