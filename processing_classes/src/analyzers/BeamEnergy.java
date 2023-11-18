
package analyzers;

/**
 *
 * @author tbhayward
 */
public class BeamEnergy {
    protected double Eb;
    
    public BeamEnergy(int runnum) {
        // default beam energy set to rga fall 2018
        Eb = 10.6041; // rga fall 2018
        if (runnum >= 5032 && runnum <= 5666) { Eb = 10.6041; } // RGA Fall 2018
        else if (runnum >= 5875 && runnum <= 5995) { Eb = 6.535; } // RGK Fall 2018
        else if (runnum >= 5674 && runnum <= 5870) { Eb = 7.546; } // RGK Fall 2018
        else if (runnum >= 6616 && runnum <= 6783) { Eb = 10.1998; } // RGA Spring 2019
        else if (runnum >= 6120 && runnum <= 6399) { Eb = 10.5986; } // RGB Spring 2019
        else if (runnum >= 6409 && runnum <= 6604) { Eb = 10.1998; } // RGB Spring 2019
        else if (runnum >= 11093 && runnum <= 11283) { Eb = 10.4096; } // RGB Fall 2019
        else if (runnum >= 11284 && runnum <= 11300) { Eb = 4.17179; } // RGB Fall 2019
        else if (runnum >= 11323 && runnum <= 11571) { Eb = 10.3894; } // RGB Spring 2020
        else if (runnum >= 16290 && runnum <= 16353) { Eb = 10.5473; } // RGC preliminary summer 2022
    }
    
    public double Eb() { return Eb;}
    
}
