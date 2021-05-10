/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package org.ualg.pt.pmfambiguity;

import java.util.ArrayList;
import org.ualg.pt.pmfcore.DigesterAndMassCalculater;

/**
 * Information of ambiguous symbol Z
 * @author Eman
 */
public class ZSymbol extends AmbiguousSymbol {

    public ZSymbol(int index) {
        super();
        clearName = "Z symbol";
        charcode = 'Z';
        averageMass = 146.6375;
        correspondingChars = new ArrayList(DigesterAndMassCalculater.getZsymbols());
        this.positionInSequence=index;
    }
}
