/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package org.ualg.pt.pmfcore;

/**
 * Get report that gives an information of MS masses, contaminants database, and query database.
 * @author Eman
 */
public  class ReportModel {
    
    int proteinsNO;
    long dBsize;
    int XNO;
    double  maxProteinMass;
    double maxPeptideMass;
    

    /**
     * Get The number of proteins in the database
     * @return
     */
    public int getProteinsNO() {
        return proteinsNO;
    }

    /**
     * Set the number of proteins in the database
     * @param proteinsNO
     */
    public void setProteinsNO(int proteinsNO) {
        this.proteinsNO = proteinsNO;
    }

    /**
     * Get the size of the database
     * @return protein database size
     */
    public long getdBsize() {
        return dBsize;
    }

    /**
     *  Set database size
     * @param dBsize
     */
    public void setdBsize(long dBsize) {
        this.dBsize = dBsize;
    }

    /**
     * Get the number of X in the database
     * @return number of Xs
     */
    public int getXNO() {
        return XNO;
    }

    /**
     * Set the number of X in the database
     * @param XNO
     */
    public void setXNO(int XNO) {
        this.XNO = XNO;
    }

    /**
     * Get the maximum protein mass in the database 
     * @return maximum protein mass
     */
    public double getMaxProteinMass() {
        return maxProteinMass;
    }

    /**
     * Set the maximum protein mass in the database 
     * @param maxProteinMass
     */
    public void setMaxProteinMass(double maxProteinMass) {
        this.maxProteinMass = maxProteinMass;
    }

    /**
     * Get the maximum mass of a peptide in database
     * @return maximum peptide mass
     */
    public double getMaxPiptideMass() {
        return maxPeptideMass;
    }

    /**
     * Set the maximum mass of a peptide in database
     * @param maxPeptideMass
     */
    public void setMaxPiptideMass(double maxPeptideMass) {
        this.maxPeptideMass = maxPeptideMass;
    }   

}
