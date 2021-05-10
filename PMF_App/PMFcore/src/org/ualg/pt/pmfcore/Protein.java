/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package org.ualg.pt.pmfcore;

import java.io.File;
import java.util.ArrayList;
import java.util.Iterator;
import java.util.LinkedList;
import java.util.List;
import org.biojava.bio.BioException;
import org.biojavax.bio.seq.RichSequence;

/**
 *
 * @author Eman
 */
public class Protein {
    ArrayList<Double> massesList;
    File file;
    List ls;
    public enum ProteinType{EXPERIMENTAL,THEORETICAL,CONTAMINANTS};
    ProteinType type;
    String ProteinID;
    String Name;
    RichSequence sequence;
    double proteinMass;
    double proteinPI;
    int numberOfMatched;
    LinkedList<Peptide> matchedPeptides=new LinkedList<Peptide>();
    LinkedList<Peptide> peptides=new LinkedList<Peptide>();
    
    public Protein(Protein p) throws BioException {
       this.type=p.getType();
       if(p.getMassesList()!=null){
       this.massesList= new ArrayList<Double>();
       for(Iterator it=p.getMassesList().iterator();it.hasNext();) {
               this.massesList.add((Double) it.next());
           }
       }
       
       if(p.getSequence()!=null){
           
           this.sequence= RichSequence.Tools.createRichSequence(p.getSequence().getNamespace(),p.getSequence().getName(),p.getSequence().seqString(),p.getSequence().getAlphabet());
       }
         
       this.ProteinID=p.getProteinID();
       this.Name=p.Name;
       this.proteinPI=p.proteinPI;
       this.file=p.file;
       
    }

    public Protein() {
    }

    
    /**
     * Get protein isoelectric point (PI)
     * @return protein 
     */
    public double getProteinPI() {
        return proteinPI;
    }

    /**
     * Set protein isoelectric point (PI)
     * @param proteinPI
     */
    public void setProteinPI(double proteinPI) {
        this.proteinPI = proteinPI;
    }
    

    /**
     * Get the number of matched peptides between experimental sample and database protein.
     * @return number of matches
     */
    public int getNumberOfMatched() {
        return numberOfMatched;
    }

    /**
     * Set the number of matched peptides between experimental sample and database protein.
     * @param numberOfMatched
     */
    public void setNumberOfMatched(int numberOfMatched) {
        this.numberOfMatched = numberOfMatched;
    }

    
    /**
     * Get the peptide sequence that its mass matches the peak list mass
     * @return peptide sequence
     */
    public LinkedList<Peptide> getMatchedPeptides() {
        return matchedPeptides;
    }

    /**
     * Set the peptide sequence that its mass matches the peak list mass
     * @param matchedPeptides
     */
    public void setMatchedPeptides(LinkedList<Peptide> matchedPeptides) {
        this.matchedPeptides = matchedPeptides;
    }

    /**
     * Get peptides sequences of protein
     * @return peptides sequences
     */
    public LinkedList<Peptide> getPeptides() {
        return peptides;
    }

    /**
     * Set peptides sequences 
     * @param unmatchedPeptides
     */
    public void setPeptides(LinkedList<Peptide> unmatchedPeptides) {
        this.peptides = unmatchedPeptides;
    }
    
    
    /**
     * Get protein name
     * @return protein name
     */
    public String getName() {
        return Name;
    }

    /**
     * set protein name
     * @param Name
     */
    public void setName(String Name) {
        this.Name = Name;
    }
    
   
    /**
     * Get protein molecular weight
     * @return protein mass
     */
    public double getProteinMass() {
        return proteinMass;
    }

    /**
     * set protein molecular weight
     * @param proteinMass
     */
    public void setProteinMass(double proteinMass) {
        this.proteinMass = proteinMass;
    }

    /**
     * Get protein type: theoretical or experimental
     * @return protein type
     */
    public ProteinType getType() {
        return type;
    }

    /**
     * Set protein type: theoretical or experimental 
     * @param type
     */
    public void setType(ProteinType type) {
        this.type = type;
    }
                  
    
    /**
     * get a list of peptides masses of  digested protein 
     * @return masses list
     */
    public ArrayList<Double> getMassesList() {
        return massesList;
    }

    /**
     * Set peptides masses of digested protein
     * @param massesList
     */
    public void setMassesList(ArrayList<Double> massesList) {
        this.massesList = massesList;
    }

    /**
     * Get file where the protein is existed 
     * @return protein file 
     */
    public File getFile() {
        return file;
    }

    /**
     * Set protein file
     * @param file
     */
    public void setFile(File file) {
        this.file = file;
    }

    /**
     * get protein identifier 
     * @return protein ID
     */
    public String getProteinID() {
        return ProteinID;
    }

    /**
     * Set protein identifier
     * @param ProteinID
     */
    public void setProteinID(String ProteinID) {
        this.ProteinID = ProteinID;
    }

    /**
     * Get protein sequence
     * @return protein sequence
     */
    public RichSequence getSequence() {
        return sequence;
    }

    /**
     * Set protein sequence
     * @param seq
     */
    public void setSequence(RichSequence seq) {
        this.sequence = seq;
    }
    
    /**
     * Check if the protein is digested
     * @return
     */
    public boolean isDigested(){
        
        return massesList!=null;
    }

   /**
     * Get protein information
     * @return protein information
     */ 
    @Override
    public String toString() {

        String massesString="[";
        StringBuilder proteinSequence= new StringBuilder(massesString);
        
        if(massesList!=null){
            Iterator it= massesList.iterator();
        
        while(it.hasNext()){
            double massValue= (Double)it.next();
           // massesString+= String.valueOf(massValue);
            proteinSequence.append(massValue);
            if(it.hasNext()) {
                //massesString+=",";
                proteinSequence.append(",");
            }
        }
        }
        massesString+="]";
            
        String peptides="";
        for(Peptide p: matchedPeptides) {
            peptides+=p;
        }
        
        return "Protein{" + " file=" + file + ", type=" + type + ", ProteinID=" + ProteinID + ", sequence=" + sequence.seqString() + ", massesList=" + proteinSequence+']'+ "\n Matched Matches:\n"+ peptides+  '}';
    } 
    
    
    /**
     * Get matched peptides information: sequence, START position, and END position
     * @return peptide
     */
    public String matchedPeptideAsString(){
       String peptides="";
       double adjustedMass;
        for(Peptide p: matchedPeptides){
            adjustedMass=(double)Math.round(p.getMass()*1000)/1000;
            peptides+= "<font color=\"red\">"+ p.peptideValue + "</font>"+" &nbsp \t" + p.getStart() + " ... &nbsp" + p.getEnd() + "&nbsp \t"+adjustedMass + "\n <br />"    ;
        }
        return peptides;
   } 
   
    public String proteinSequenceWithMatched(Protein protein){
        String sequence = protein.getSequence().seqString();
        String newString = "<table> <tr><td class=\"from\">";
            for(int i=0; i<protein.peptides.size();i++){
            for(Peptide matchedPeptide: matchedPeptides){
                if(protein.peptides.get(i).getStart()==matchedPeptide.getStart())
                protein.peptides.get(i).peptideValue= "<font color=\"red\">"+ protein.peptides.get(i).peptideValue+"</font>";              
        }
          
            
            }  
            for(int i=0; i<protein.peptides.size();i++)
                 newString+=protein.peptides.get(i).peptideValue;
            newString+="</td> </tr> </table>";
        return newString;
    }
    // statistics
    
    /**
     * Get the matched peptides length coverage
     * @return length coverage
     */
    public double getCoverage(){
        double length=0;
        for (Peptide p: matchedPeptides){
            length+=p.getPepTideValue().length();
        }
        return length/(double)sequence.seqString().length();
    }
    
    
    /**
     * Get the matched coverage ratio
     * @return matched coverage
     */
    public double getMatchedCoverage(){
        double correctedNumber=(double) getCorrectedNumber(2);
        if(correctedNumber==0) {
            return 0;
        }
        return (double )matchedPeptides.size()/ correctedNumber;
    }
    
    /**
     * ignore the small peptides that only contain one or two amino acids 
     * @param ignorlength
     * @return amino acids number in the peptide to be considered
     */
    private int getCorrectedNumber(int ignorlength){
        int number=0;
        for(Peptide p: peptides) {
            if(p.getPeptideLength()>ignorlength) {
                number++;
            }
        }
            return number;
    }
    
    /**
     * Get the matched peptide mass coverage
     * @return  mass coverage
     */
    public double getCoverageMass(){
        double masses=0;
        for (Peptide p: matchedPeptides){
            masses+=p.getMass();
        }
        return masses/(double)this.proteinMass;
    }
    
    /**
     * Add matched peptides to the matchedPeptides list without duplicate. 
     * @param p
     */
    public void addMatchedPeptide(Peptide p){
        for(Peptide peptide: matchedPeptides){
            if(p.getPepTideValue().equals(peptide.getPepTideValue())) {
                return;
            } // duplicated reject
            
        }
        matchedPeptides.addLast(p);
    } 
//    public StringBuilder getProteinSequence(){
//        StringBuilder proteinString= new StringBuilder();
//        
//        proteinString.append("<table> <tr><td>");
//        Iterator peptideIt=sequence.iterator();
//        Peptide onePeptide= (Peptide)peptideIt.next();
//        
//        while(peptideIt.hasNext()){
//            String matchedPeptide=matchedPeptideAsString();
//            if(onePeptide.getPepTideValue().equals(matchedPeptide))
//                proteinString.append("<font color=\"red\">")
//                        .append(onePeptide).append("</font> " );
//                
//        }
//        proteinString.append("</td> </tr> </table>");
//        
//        return proteinString;
//        
//    }
}
