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
public class MatchComparator implements Comparator<Match> {
    
      public int compare(Match m1, Match m2){
            if (m1.getScore() < m2.getScore()) return 1;
            if (m1.getScore() > m2.getScore()) return -1;
            return 0;
        }
    
}
