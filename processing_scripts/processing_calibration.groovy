/*
 * author Timothy B. Hayward
 * 
 * calibration 
 */

// import CLAS12 physics classes
import org.jlab.io.hipo.*
import org.jlab.io.base.DataEvent
import org.jlab.clas.physics.*
import org.jlab.clas12.physics.*

// filetype for gathering files in directory
import groovy.io.FileType

// dilks CLAS QA analysis
import clasqa.QADB

class CalibrationScript {
    // Define instance variables with default values
    int config_run = -999
    int config_event = -999
    int config_trigger = -999
    double torus = -999
    double solenoid = -999

    int event_helicity = 0

    int particle_pid = -999
    double particle_px = -999
    double particle_py = -999
    double particle_pz = -999
    double p = -999;
    double particle_vx = -999
    double particle_vy = -999
    double particle_vz = -999
    double particle_beta = -999
    double particle_chi2pid = -999
    int particle_status = -999
    double theta = -999
    double phi = -999

    int cal_sector = -999
    double cal_energy_1 = -999; double cal_energy_4 = -999; double cal_energy_7 = -999;
    double cal_x_1 = -999; double cal_x_4 = -999; double cal_x_7 = -999; 
    double cal_y_1 = -999; double cal_y_4 = -999; double cal_y_7 = -999;
    double cal_z_1 = -999; double cal_z_4 = -999; double cal_z_7 = -999;
    double cal_lu_1 = -999; double cal_lu_4 = -999; double cal_lu_7 = -999;
    double cal_lv_1 = -999; double cal_lv_4 = -999; double cal_lv_7 = -999;
    double cal_lw_1 = -999; double cal_lw_4 = -999; double cal_lw_7 = -999;

    int cc_sector = -999
    double cc_nphe_15 = -999; double cc_nphe_16 = -999; 

    int track_sector_5 = -999; int track_sector_6 = -999; 
    double track_chi2_5 = -999; double track_chi2_6 = -999; 
    int track_ndf_5 = -999; int track_ndf_6 = -999; 

    int traj_detector = -999
    int traj_layer = -999
    double traj_x = -999
    double traj_y = -999
    double traj_z = -999
    double traj_edge = -999

    // Method to reset all variables to their default values
    void resetVariables() {
        // config_run = -999
        // config_event = -999
        // config_trigger = -999
        // torus = -999
        // solenoid = -999

        // event_helicity = 0

        particle_pid = -999;
	    particle_px = -999;
	    particle_py = -999;
	    particle_pz = -999;
	    p = -999;
	    particle_vx = -999;
	    particle_vy = -999;
	    particle_vz = -999;
	    particle_beta = -999;
	    particle_chi2pid = -999;
	    particle_status = -999;
	    theta = -999;
	    phi = -999;

	    cal_sector = -999;
	    cal_energy_1 = -999; cal_energy_4 = -999; cal_energy_7 = -999;
	    cal_x_1 = -999; cal_x_4 = -999; cal_x_7 = -999; 
	    cal_y_1 = -999; cal_y_4 = -999; cal_y_7 = -999;
	    cal_z_1 = -999; cal_z_4 = -999; cal_z_7 = -999;
	    cal_lu_1 = -999; cal_lu_4 = -999; cal_lu_7 = -999;
	    cal_lv_1 = -999; cal_lv_4 = -999; cal_lv_7 = -999;
	    cal_lw_1 = -999; cal_lw_4 = -999; cal_lw_7 = -999;

	    cc_sector = -999;
	    cc_nphe_15 = -999; cc_nphe_16 = -999; 

	    track_sector_5 = -999; track_sector_6 = -999; 
	    track_chi2_5 = -999; track_chi2_6 = -999; 
	    track_ndf_5 = -999; track_ndf_6 = -999; 

	    traj_detector = -999;
	    traj_layer = -999;
	    traj_x = -999;
	    traj_y = -999;
	    traj_z = -999;
	    traj_edge = -999;
    }

    // Helper method for formatting doubles
    String formatDouble(double value) {
    	if (value == -999.0) { return -999; }
        else { return String.format("%.3f", value) }
    }

    // Static method to calculate phi
	static double phi_calculation(double x, double y) {
	    double phi = Math.toDegrees(Math.atan2(x, y))
	    phi = phi - 90
	    if (phi < 0) {
	        phi = 360 + phi
	    }
	    phi = 360 - phi
	    return phi
	}

	// Static method to calculate theta
	static double theta_calculation(double x, double y, double z) {
	    double r = Math.sqrt(x * x + y * y + z * z)
	    return Math.toDegrees(Math.acos(z / r))
	}

    // Method for the main logic
    void run(String[] args) {
        // Start time
        long startTime = System.currentTimeMillis()

        // ~~~~~~~~~~~~~~~~ set up input parameters ~~~~~~~~~~~~~~~~ //

        // Check if an argument is provided
        if (!args) {
            // Print an error message and exit the program if the input directory is not specified
            println("ERROR: Please enter a hipo file directory as the first argument")
            System.exit(0)
        }
        // If the input directory is provided, iterate through each file recursively
        def hipo_list = []
        (args[0] as File).eachFileRecurse(FileType.FILES) { if (it.name.endsWith('.hipo')) hipo_list << it }

        // Set the output file name based on the provided 3rd argument or use the default name
        String output_file = args.length < 2 ? "inclusive_dummy_out.txt" : args[1]
        if (args.length < 2)
            println("WARNING: Specify an output file name. Set to \"hadron_dummy_out.txt\".")
        File file = new File(output_file)
        file.delete()
        BufferedWriter writer = new BufferedWriter(new FileWriter(file))

        // Set the number of files to process based on the provided 4th argument
        // use the size of the hipo_list if no argument provided
        int n_files = args.length < 3 || Integer.parseInt(args[2]) > hipo_list.size()
            ? hipo_list.size() : Integer.parseInt(args[2])
        if (args.length < 3 || Integer.parseInt(args[2]) > hipo_list.size()) {
            // Print warnings and information if the number of files is not specified or too large
            println("WARNING: Number of files not specified or number too large.")
            println("Setting # of files to be equal to number of files in the directory.")
            println("There are $hipo_list.size files.")
        }

        // ~~~~~~~~~~~~~~~~ prepare physics analysis ~~~~~~~~~~~~~~~~ //

        // setup QA database
        QADB qa = new QADB();

        // create a StringBuilder for accumulating lines
        StringBuilder batchLines = new StringBuilder()

        int num_events = 0
        int max_lines = 1000
        int lineCount = 0
        for (current_file in 0..<n_files) {
            // limit to a certain number of files defined by n_files
            println("\n Opening file "+Integer.toString(current_file+1)
                +" out of "+n_files+".\n")

            HipoDataSource reader = new HipoDataSource()
            reader.open(hipo_list[current_file]) // open next hipo file
            HipoDataEvent event = reader.getNextEvent()

            while (reader.hasEvent()) {
                ++num_events
                if (num_events % 500000 == 0) { // not necessary, just updates output
                    print("processed: " + num_events + " events. ")
                }
                // get run and event numbers
                event = reader.getNextEvent()

                HipoDataBank run_Bank = (HipoDataBank) event.getBank("RUN::config")
                HipoDataBank event_Bank = (HipoDataBank) event.getBank("REC::Event")
                HipoDataBank rec_Bank = (HipoDataBank) event.getBank("REC::Particle")
                HipoDataBank cal_Bank = (HipoDataBank) event.getBank("REC::Calorimeter")
                HipoDataBank cc_Bank = (HipoDataBank) event.getBank("REC::Cherenkov")
                HipoDataBank track_Bank = (HipoDataBank) event.getBank("REC::Track")
                HipoDataBank traj_Bank = (HipoDataBank) event.getBank("REC::Traj")

                // collect info for QA
                config_run = run_Bank.getInt('run', 0)
                config_event = run_Bank.getInt('event', 0)

                // do not use the qa if it is MC (runnum = 11)
                // do not use the qa if the run is from RGC (until QA is produced!)
                // boolean process_event = filter.isValid(research_Event)
                boolean process_event = (config_run == 11 || config_run < 5020 ||
                    config_run >= 11571 || qa.OkForAsymmetry(config_run, config_event))

                if (process_event) {

                    event_helicity = event_Bank.getInt('helicity',0);

                    for (int particle_Index = 0; particle_Index < rec_Bank.rows(); 
                    	particle_Index++) {

                        particle_pid = rec_Bank.getInt("pid", particle_Index);
                        if (particle_pid == 0 || particle_pid == 45) { continue; }
                        particle_px = rec_Bank.getFloat("px", particle_Index);
                        particle_py = rec_Bank.getFloat("py", particle_Index);
                        particle_pz = rec_Bank.getFloat("pz", particle_Index);
						particle_vx = rec_Bank.getFloat("vx",particle_Index);
						particle_vy = rec_Bank.getFloat("vy",particle_Index);
						particle_vz = rec_Bank.getFloat("vz",particle_Index);
						p = Math.sqrt(particle_px*particle_px+
							particle_py*particle_py+particle_pz*particle_pz);
						theta = theta_calculation(particle_px, particle_py, particle_pz);
						phi = phi_calculation(particle_px, particle_py);
						particle_beta = rec_Bank.getFloat("beta",particle_Index);
						particle_chi2pid = rec_Bank.getFloat("chi2pid",particle_Index);
						particle_status = rec_Bank.getInt("status",particle_Index);

						// Calorimeter
						for (int current_Row = 0; current_Row < cal_Bank.rows(); current_Row++) {
						    int pindex = cal_Bank.getInt("pindex", current_Row);
						    if (pindex == particle_Index) {
						        cal_sector = cal_Bank.getInt("sector", current_Row);
						        int layer = cal_Bank.getInt("layer", current_Row);
						        float energy = cal_Bank.getFloat("energy", current_Row);
						        float x = cal_Bank.getFloat("x", current_Row);
						        float y = cal_Bank.getFloat("y", current_Row);
						        float z = cal_Bank.getFloat("z", current_Row);
						        float lu = cal_Bank.getFloat("lu", current_Row);
						        float lv = cal_Bank.getFloat("lv", current_Row);
						        float lw = cal_Bank.getFloat("lw", current_Row);

						        switch (layer) {
						            case 1:
						                cal_energy_1 = energy;
						                cal_x_1 = x;
						                cal_y_1 = y;
						                cal_z_1 = z;
						                cal_lu_1 = lu;
						                cal_lv_1 = lv;
						                cal_lw_1 = lw;
						                break;
						            case 4:
						                cal_energy_4 = energy;
						                cal_x_4 = x;
						                cal_y_4 = y;
						                cal_z_4 = z;
						                cal_lu_4 = lu;
						                cal_lv_4 = lv;
						                cal_lw_4 = lw;
						                break;
						            case 7:
						                cal_energy_7 = energy;
						                cal_x_7 = x;
						                cal_y_7 = y;
						                cal_z_7 = z;
						                cal_lu_7 = lu;
						                cal_lv_7 = lv;
						                cal_lw_7 = lw;
						                break;
						        }
						    }
						}


	                    // // Cherenkov Counter
	                    // for (int current_Row = 0; current_Row < cc_Bank.rows(); current_Row++) {
	                    //     // Get the pindex and layer values for the current row
	                    //     int pindex = cc_Bank.getInt("pindex", current_Row);
	                    //     if (pindex == particle_Index) {
	                    //         int cc_sector = cc_Bank.getInt("sector", current_Row);
	                    //         int detector = cc_Bank.getInt("detector", current_Row);
	                    //         double nphe = cc_Bank.getFloat("nphe", current_Row);
	                    //         switch(detector) {
	                    //         	case 15:
	                    //         		cc_nphe_15 = nphe;
	                    //         		break;
	                    //         	case 16:
	                    //         		cc_nphe_16 = nphe;
		                //             	break
	                    //         }
	                    //     }
	                    // }


	                    // // Track Bank
	                    // for (int current_Row = 0; current_Row < track_Bank.rows(); current_Row++) {
	                    //     // Get the pindex and layer values for the current row
	                    //     int pindex = track_Bank.getInt("pindex", current_Row);
	                    //     if (pindex == particle_Index) {
	                    //         int detector = track_Bank.getInt("detector", current_Row);
	                    //         int sector = track_Bank.getInt("sector", current_Row);
	                    //         double chi2 = track_Bank.getFloat("chi2", current_Row);
	                    //         int ndf = track_Bank.getInt("NDF", current_Row);
	                    //         switch(detector) {
	                    //         	case 5:
	                    //         		track_sector_5 = sector;
	                    //         		track_chi2_5 = chi2;
	                    //         		track_ndf_5 = ndf;
	                    //         		break;
	                    //         	case 6:
	                    //         		track_sector_6 = sector;
	                    //         		track_chi2_6 = chi2;
	                    //         		track_ndf_6 = ndf;
		                //             	break
	                    //         }
	                    //     }
	                    // }

	                    // // Use a StringBuilder to append all data in a single call
	                    // StringBuilder line = new StringBuilder()
	                    // line.append(config_run).append(" ")
	                    //     .append(config_event).append(" ")

	                    //     .append(event_helicity).append(" ")

	                    //     .append(particle_pid).append(" ")
	                    //     .append(formatDouble(particle_px)).append(" ")
	                    //     .append(formatDouble(particle_py)).append(" ")
	                    //     .append(formatDouble(particle_pz)).append(" ")
	                    //     .append(formatDouble(p)).append(" ")
	                    //     .append(formatDouble(theta)).append(" ")
	                    //     .append(formatDouble(phi)).append(" ")
	                    //     .append(formatDouble(particle_vx)).append(" ")
	                    //     .append(formatDouble(particle_vy)).append(" ")
	                    //     .append(formatDouble(particle_vz)).append(" ")
	                    //     .append(formatDouble(particle_beta)).append(" ")
	                    //     .append(formatDouble(particle_chi2pid)).append(" ")
	                    //     .append(particle_status).append(" ")

	                    //     .append(cal_sector).append(" ")
	                    //     .append(formatDouble(cal_energy_1)).append(" ")
	                    //     .append(formatDouble(cal_x_1)).append(" ")
	                    //     .append(formatDouble(cal_y_1)).append(" ")
	                    //     .append(formatDouble(cal_z_1)).append(" ")
	                    //     .append(formatDouble(cal_lu_1)).append(" ")
	                    //     .append(formatDouble(cal_lv_1)).append(" ")
	                    //     .append(formatDouble(cal_lw_1)).append(" ")
	                    //     .append(formatDouble(cal_x_4)).append(" ")
	                    //     .append(formatDouble(cal_y_4)).append(" ")
	                    //     .append(formatDouble(cal_z_4)).append(" ")
	                    //     .append(formatDouble(cal_lu_4)).append(" ")
	                    //     .append(formatDouble(cal_lv_4)).append(" ")
	                    //     .append(formatDouble(cal_lw_4)).append(" ")
	                    //     .append(formatDouble(cal_x_7)).append(" ")
	                    //     .append(formatDouble(cal_y_7)).append(" ")
	                    //     .append(formatDouble(cal_z_7)).append(" ")
	                    //     .append(formatDouble(cal_lu_7)).append(" ")
	                    //     .append(formatDouble(cal_lv_7)).append(" ")
	                    //     .append(formatDouble(cal_lw_7)).append(" ")

	                    //     .append(cc_sector).append(" ")
	                    //     .append(formatDouble(cc_nphe_15)).append(" ")
	                    //     .append(formatDouble(cc_nphe_16)).append(" ")

	                    //     .append(track_sector_5).append(" ")
	                    //     .append(formatDouble(track_chi2_5)).append(" ")
	                    //     .append(track_ndf_5).append(" ")
	                    //     .append(track_sector_6).append(" ")
	                    //     .append(formatDouble(track_chi2_6)).append(" ")
	                    //     .append(track_ndf_6).append(" ")

	                    //     .append(config_run).append("\n")

	                    // // Append the line to the batchLines StringBuilder
	                    // batchLines.append(line.toString())
	                    // lineCount++ // Increment the line count

	                    // // If the line count reaches 1000, write to the file and reset
	                    // if (lineCount >= max_lines) {
	                    //     file.append(batchLines.toString())
	                    //     batchLines.setLength(0)
	                    //     lineCount = 0
	                    // }

	                    // // Reset variables after processing each particle
	                    // resetVariables()
	                }
	            }
            	reader.close()
	        }

	        // Write any remaining lines in the batchLines StringBuilder to the file
	        if (batchLines.length() > 0) {
	            file.append(batchLines.toString())
	            batchLines.setLength(0)
	        }

	        // println("\n1:runnum, 2:evnum, 3:helicity, 4:e_p, 5:e_theta, 6:e_phi, 7:vz_e,"+
	        // "8:Q2, 9:W, 10:Mx, 11: Mx2, 12:x, 13:y, 14: DepA, 15: DepB, 16: DepC, 17: DepV, 18: DepW\n")

	        println("output text file is: $file")
	    }

	    writer.close()

	    // End time
	    long endTime = System.currentTimeMillis()
	    // Calculate the elapsed time
	    long elapsedTime = endTime - startTime
	    // Print the elapsed time in milliseconds
	    println("Elapsed time: ${elapsedTime} ms")
	}
}

// Create an instance of the script and run it
def script = new CalibrationScript()
script.run(args)