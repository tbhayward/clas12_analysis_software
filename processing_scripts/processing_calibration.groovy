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

// def banks_test = { DataEvent event ->
//     String[] bankNames = 
//         ["RUN::config","REC::Event","REC::Particle","REC::Calorimeter",
//          "REC::Track","REC::Traj","REC::Cherenkov"]
//     for (String bankName : bankNames) {
//         if (!event.hasBank(bankName)) { return false }
//     }
//     return true
// }

public static void main(String[] args) {

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

    // RUN::config variables
    int config_run = -9999; int config_event = -9999; int config_trigger = -9999
    double torus = -9999; double solenoid = -9999
    // double is correct, a few runs in RGA are torus 1.06 instead of +1

    // REC::Event variables
    int event_helicity = 0

    // REC::Particle variables
    int particle_pid = -9999
    double particle_px = -9999; double particle_py = -9999; double particle_pz = -9999
    double particle_vx = -9999; double particle_vy = -9999; double particle_vz = -9999
    double particle_beta = -9999; double particle_chi2pid = -9999
    int particle_status = -9999

    // REC::Calorimeter variables
    int cal_sector = -9999
    int cal_layer = -9999
    double cal_energy = -9999
    double cal_x = -9999; double cal_y = -9999; double cal_z = -9999
    double cal_lu = -9999; double cal_lv = -9999; double cal_lw = -9999

    // REC::Cherenkov variables
    int cc_sector = -9999
    int cc_layer = -9999
    double cc_nphe = -9999

    // REC::Track variables
    int track_detector = -9999; int track_sector = -9999
    double track_chi2 = -9999
    int track_ndf = -9999

    // REC::Traj variables
    int traj_detector = -9999; int traj_layer = -9999
    double traj_x = -9999; double traj_y = -9999; double traj_z = -9999
    double traj_edge = -9999

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
                
                // Use a StringBuilder to append all data in a single call
                StringBuilder line = new StringBuilder()
                line.append(config_run).append(" ")
                    .append(config_event).append(" ")
                    .append(Depolarization_W).append("\n")

                // Append the line to the batchLines StringBuilder
                batchLines.append(line.toString())
                lineCount++ // Increment the line count

                // If the line count reaches 1000, write to the file and reset
                if (lineCount >= max_lines) {
                    file.append(batchLines.toString())
                    batchLines.setLength(0)
                    lineCount = 0
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