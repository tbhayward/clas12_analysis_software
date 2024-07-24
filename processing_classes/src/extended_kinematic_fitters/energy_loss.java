/**
 *
 * @author Timothy B. Hayward
 */

package extended_kinematic_fitters;

import org.jlab.clas.physics.GenericKinematicFitter;
import org.jlab.clas.physics.Particle;
import org.jlab.clas.physics.PhysicsEvent;
import org.jlab.io.base.DataEvent;
import org.jlab.io.hipo.HipoDataBank;

import org.jlab.clas.physics.*;

/**
 *
 * @author tbhayward
 */
public class energy_loss {
    
    public double pass2_fd_energy_loss(float polarity, int pid, double px, double py, double pz) {
        double dp = 0; // scale size
        double r = Math.pow(Math.pow(px,2)+Math.pow(py,2)+Math.pow(pz,2),0.5);
        double theta = (180/Math.PI)*Math.acos(pz/r);
        double p = Math.sqrt(px*px+py*py+pz*pz);
        
        if (polarity < 0) {
            if (theta < 27) {
                switch (pid) { 
                    case 11:
                        dp = 0.00730953+Math.exp(-14.0303+0.957154*p);
                        break;
                        
                    case 211:
                        dp = -36.3456+Math.exp(3.59313-0.0000065*p);
                        break;

                    case -211:
                        dp = 0.00396596+Math.exp(6.14785-16.8218*p);
                        break;

                    case 321:
                        dp = -34.5647+Math.exp(3.54297-0.0000136581*p);
                        break;

                    case -321:
                        dp = 0.00340915+Math.exp(-0.567643-6.31624*p);
                        break;
                        
                    case 2212:
                        dp = 0.00428181+Math.exp(-2.62675-3.39989*p);
                        break;
                }
            } else if (theta >= 27) {
                switch (pid) { 
                    case 11:
                        dp = 0.00552679+Math.exp(-4.2381-0.660486*p);
                        break;
                        
                    case 211:
                        dp = -45.9594+Math.exp(3.82795-0.0000358352*p);
                        break;

                    case -211:
                        dp = -0.000914991+Math.exp(-4.29887-0.555948*p);
                        break;

                    case 321:
                        dp = 0.00514143+Math.exp(-3.65383-2.19434*p);
                        break;

                    case -321:
                        dp = 0.00106601+Math.exp(-3.87516-1.04774*p);
                        break;
                        
                    case 2212:
                        dp = 0.00597654+Math.exp(-2.4691-2.37887*p);
                        break;
                }
            }
        } else if (polarity > 0) {
            if (theta < 27) {
                switch (pid) { 
                    case 11:
                        dp = -30.2+Math.exp(3.40808-0.0000048*p);
                        break;
                    
                    case 211:
                        dp = 0.00392301+Math.exp(3.83445-13.9503*p);
                        break;

                    case -211:
                        dp = -39.0688+Math.exp(3.66543-0.000006173*p);
                        break;

                    case 321:
                        dp = 0.00377599+Math.exp(-0.386989-7.33292*p);
                        break;

                    case -321:
                        dp = 0.00377308+Math.exp(-2.11549-6.31415*p);
                        break;
                        
                    case 2212:
                        dp = 0.00443825+Math.exp(-0.883026-5.07364*p);
                        break;
                }
            } else if (theta >= 27) {
                switch (pid) { 
                    case 11:
                        dp = -0.0154673+Math.exp(-3.42119-0.096196*p);
                        break;
                        
                    case 211:
                        dp = 0.00185544+Math.exp(-4.41637-0.878962*p);
                        break;

                    case -211:
                        dp = -91.7796+Math.exp(4.51949-0.0000206923*p);
                        break;

                    case 321:
                        dp = 0.00279201+Math.exp(-3.64049-1.65387*p);
                        break;

                    case -321:
                        dp = 0.00430273+Math.exp(-3.86056-1.6304*p);
                        break;
                        
                    case 2212:
                        dp = 0.00475655+Math.exp(-1.99834-2.81217*p);
                        break;
                }
            }
        }
        
        double fe = (dp+p)/p;
        return fe;
        
    }
}
