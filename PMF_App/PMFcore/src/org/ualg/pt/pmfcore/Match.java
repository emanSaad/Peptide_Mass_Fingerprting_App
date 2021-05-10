/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package org.ualg.pt.pmfcore;

import java.util.Comparator;

/**
 *
 * @author faroq
 */
    public class Match{

        String key;
        double score;
        Protein protein;
        double expProteins;
        double maxScore;
        double minScore;
        double significance;

        public String getKey() {
            return key;
        }
        public Match(String key, double score, Protein protein, double exProteins) {
            this.key = key;
            this.score = score;
            this.protein = protein;
            this.expProteins = exProteins;
        }

        public void setKey(String key) {
            this.key = key;
        }
        
        public double getMaxScore( ){ 
            return maxScore;
            
        }
        
        public void setMaxScore(double MaxScore){
            this.maxScore=maxScore; 
        }
        public double getMinScore( ){ 
            return minScore;
            
        }
        
        public void setMinScore(double MinScore){
            this.minScore=minScore; 
        }
        

        public Match copy() {
            return new Match(key, score, protein, expProteins);
        }

        public double getScore() {
            return score;
        }

        public void setScore(double score) {
            this.score = score;
        }

        public Protein getProtein() {
            return protein;
        }

        public void setProtein(Protein protein) {
            this.protein = protein;
        }

        @Override
        public String toString() {
            
            //DecimalFormat df = new DecimalFormat("#");
             //df.setMaximumFractionDigits(0);
            //DecimalFormat eFormat = new DecimalFormat("0.##E0");
         
            //String proScore = eFormat.format(score);
           //String simpleScore=df.format(score);
            //double simpleScore=(double) Math.round(-Math.log10(1 / score) * 1000) / 1000;
           //double finalScore=score*10;
            return "Match{" + "Protein=" +key + "&nbsp &nbsp \t\tAccession&nbsp" + protein.getSequence().getAccession() + " &nbsp &nbsp \t\tscore=" + score +" &nbsp &nbsp \t\t Matched=" + protein.getMatchedPeptides().size()+ "&nbsp &nbsp \t\t score=" + score+ "}\n\n, <h4>Protein Details:</h4>\n" + "Mass:" + (double) Math.round(protein.proteinMass * 1000) / 1000 + " <br />\nLength:" + protein.getSequence().seqString().length() + "<br /> \nLength Coverage:" + (double) Math.round(protein.getCoverage() * 1000) / 1000 + "\n <br /> Matched Ratio:" + (double) Math.round(protein.getMatchedCoverage() * 1000) / 1000 + /*"\nSignificance:"+ (double) Math.round((protein.getCoverage()* protein.getMatchedCoverage()* ((double)protein.numberOfMatched/expProteins))* 1000)/1000 +*/ "\n <br /> Matched Masses Coverage: " + (double) Math.round(protein.getCoverageMass() * 1000) / 1000 + "\n <br /> Matched Peptides: <br />\n" + protein.matchedPeptideAsString() + "<br /> Protein Sequence: "+ protein.proteinSequenceWithMatched(protein) + " <br /><br />\n\n";
            // return "Match{" + "Protein=" + key + "\tAccession:"+ protein.getSequence().getAccession()+ "\tscore=" + ProScore  + "\t Matched=" + protein.getMatchedPeptides().size()+" \t -log (score)="+ (double) Math.round(-Math.log10(1/score) *1000)/1000 +"\n";
            //return "Match{" + "Protein=" + key + "\tAccession:"+ protein.getSequence().getAccession()+ "\tscore=" + score + "\t Matched=" + protein.getMatchedPeptides().size()+" \t -log (score)="+ (-Math.log10(1/score)) +"}\n, Protein Details:\n" + "Mass:"+ protein.proteinMass +"\nLength:" +protein.getSequence().seqString().length() + "\nLength Coverage:"+ protein.getCoverage() + "\nMatched Ratio:" + protein.getMatchedCoverage() + "\n Significance:"+ ((protein.getCoverage()* protein.getMatchedCoverage()* ((double)protein.numberOfMatched/expProteins))) +   "\nMatched Masses Coverage:"+ protein.getCoverageMass() +"\nMatched Peptides:\n"  + protein.matchedPeptideAsString()+ "\n\n";
        }       
        public String printProteinINFO(){
            //DecimalFormat df = new DecimalFormat("#");
             //df.setMaximumFractionDigits(0);
             
               
           // DecimalFormat eFormat = new DecimalFormat("0.##E0");
            //String proScore = eFormat.format(score);
            //String simpleScore=df.format(score);
            
           // double simpleScore=(double) Math.round(-Math.log10(1 / score) * 1000) / 1000;
           //double finalScore=simpleScore*10;
            
          // return   protein.getSequence().getAccession()  +"\t\t"+ proScore.toString()+"\t\t"+ /*finalScore+"\t\t"+*/ protein.getMatchedPeptides().size(); //+"</td></tr></table>"; 
          return   key  +"\t\t"+ score+"\t\t"+ /*finalScore+"\t\t"+*/ protein.getMatchedPeptides().size(); //+"</td></tr></table>"; 

        }
    }

