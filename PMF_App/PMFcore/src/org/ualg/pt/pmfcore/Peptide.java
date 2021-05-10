/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package org.ualg.pt.pmfcore;

/**
 *
 * @author Eman
 */
public class Peptide {
    String peptideValue;
    double mass;
    int start;
    int end;
    int peptideLength;
    

    /**
     * Get the peptide length
     * @return peptide length
     */
    public int getPeptideLength() {
        return end-start;
    }

    /**
     * Set the peptide length
     */
    public void setPeptideLength() {
        this.peptideLength = end-start;
    }

    /**
     * Get the sequence of amino acids for the peptide
     * @return 
     */
    public String getPepTideValue() {
        return peptideValue;
    }

    /**
     * Set the sequence of amino acids for the peptide
     * @param pepTideValue
     */
    public void setPepTideValue(String pepTideValue) {
        this.peptideValue = pepTideValue;
    }

    /**
     * Get peptide molecular weight
     * @ return peptide molecular weight
     */
    public double getMass() {
        return mass;
    }

    /**
     * Set peptide molecular weight
     * @param mass
     */
    public void setMass(double mass) {
        this.mass = mass;
    }

    /**
     * Get the start position of peptide in its sequence 
     * @return start position 
     */
    public int getStart() {
        return start;
    }

    /**
     *  Set the start position of peptide in its sequence 
     * @param start
     */
    public void setStart(int start) {
        this.start = start;
    }

    /**
     * Get the end position of peptide in its sequence 
     * @return  end position
     */
    public int getEnd() {
        return end;
    }

    /**
     * Set the end position of peptide in its sequence
     * @param end
     */
    public void setEnd(int end) {
        this.end = end;
    }

    @Override
    public String toString() {
        return "Peptide{" + "peptideValue=" + peptideValue + ", mass=" + mass + ", start=" + start + ", end=" + end + '}';
    }
    
}
