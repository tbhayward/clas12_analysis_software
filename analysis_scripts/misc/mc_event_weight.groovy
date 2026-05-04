#!/usr/bin/env groovy

import org.jlab.jnp.hipo4.io.HipoReader
import org.jlab.jnp.hipo4.data.Event
import org.jlab.jnp.hipo4.data.Bank
import org.jlab.jnp.hipo4.data.SchemaFactory

if (args.length < 1 || args.length > 2) {
    System.err.println("Usage:")
    System.err.println("  groovy sum_mc_event_weights.groovy input.hipo")
    System.err.println("  groovy sum_mc_event_weights.groovy input_directory")
    System.err.println("  groovy sum_mc_event_weights.groovy input_directory max_files")
    System.exit(1)
}

String inputPath = args[0]
int maxFiles = -1

if (args.length == 2) {
    try {
        maxFiles = Integer.parseInt(args[1])
    } catch (Exception e) {
        System.err.println("ERROR: max_files must be an integer.")
        System.exit(2)
    }

    if (maxFiles < 0) {
        System.err.println("ERROR: max_files must be non-negative.")
        System.exit(3)
    }
}

File input = new File(inputPath)

if (!input.exists()) {
    System.err.println("ERROR: input path does not exist: ${inputPath}")
    System.exit(4)
}

boolean isDirectoryInput = input.isDirectory()

List<File> hipoFiles = []

if (input.isFile()) {
    if (!input.getName().endsWith(".hipo")) {
        System.err.println("ERROR: input file does not end with .hipo: ${inputPath}")
        System.exit(5)
    }

    hipoFiles.add(input)
} else if (input.isDirectory()) {
    hipoFiles = input.listFiles()
        .findAll { File f -> f.isFile() && f.getName().endsWith(".hipo") }
        .sort { File a, File b -> a.getName() <=> b.getName() }

    if (maxFiles >= 0 && hipoFiles.size() > maxFiles) {
        hipoFiles = hipoFiles.subList(0, maxFiles)
    }
} else {
    System.err.println("ERROR: input path is neither a regular file nor a directory: ${inputPath}")
    System.exit(6)
}

if (hipoFiles.size() == 0) {
    System.err.println("ERROR: no .hipo files found to process.")
    System.exit(7)
}

long totalEvents = 0L
long totalEventsWithMcEvent = 0L
long totalRows = 0L
double totalSumWeights = 0.0
double totalMaxWeight = -Double.MAX_VALUE
boolean foundAnyWeight = false

int filesProcessed = 0

for (File hipoFile : hipoFiles) {
    println("Processing file ${filesProcessed + 1}/${hipoFiles.size()}: ${hipoFile.getPath()}")

    HipoReader reader = new HipoReader()
    reader.open(hipoFile.getPath())

    SchemaFactory factory = reader.getSchemaFactory()

    if (!factory.hasSchema("MC::Event")) {
        System.err.println("WARNING: file does not contain bank MC::Event, skipping: ${hipoFile.getPath()}")
        reader.close()
        continue
    }

    Bank mcEvent = new Bank(factory.getSchema("MC::Event"))
    Event event = new Event()

    long fileEvents = 0L
    long fileEventsWithMcEvent = 0L
    long fileRows = 0L
    double fileSumWeights = 0.0
    double fileMaxWeight = -Double.MAX_VALUE
    boolean foundFileWeight = false

    while (reader.hasNext()) {
        reader.nextEvent(event)
        event.read(mcEvent)

        fileEvents++

        int rows = mcEvent.getRows()
        if (rows <= 0) {
            continue
        }

        fileEventsWithMcEvent++

        for (int row = 0; row < rows; row++) {
            double weight = (double) mcEvent.getFloat("weight", row)

            fileSumWeights += weight
            fileRows++

            if (!foundFileWeight || weight > fileMaxWeight) {
                fileMaxWeight = weight
                foundFileWeight = true
            }

            if (!foundAnyWeight || weight > totalMaxWeight) {
                totalMaxWeight = weight
                foundAnyWeight = true
            }
        }
    }

    reader.close()

    totalEvents += fileEvents
    totalEventsWithMcEvent += fileEventsWithMcEvent
    totalRows += fileRows
    totalSumWeights += fileSumWeights
    filesProcessed++

    println("  Events read: ${fileEvents}")
    println("  Events with MC::Event: ${fileEventsWithMcEvent}")
    println("  MC::Event rows summed: ${fileRows}")
    println("  Sum of MC::Event weight: ${fileSumWeights}")

    if (foundFileWeight) {
        println("  Maximum MC::Event weight: ${fileMaxWeight}")
    } else {
        println("  Maximum MC::Event weight: none found")
    }
}

if (isDirectoryInput) {
    println("")
    println("Summary over all files looked over")
    println("Input directory: ${inputPath}")
    println("HIPO files selected: ${hipoFiles.size()}")
    println("HIPO files processed: ${filesProcessed}")
    println("Total events read: ${totalEvents}")
    println("Total events with MC::Event: ${totalEventsWithMcEvent}")
    println("Total MC::Event rows summed: ${totalRows}")
    println("Total sum of MC::Event weight: ${totalSumWeights}")

    if (foundAnyWeight) {
        println("Maximum MC::Event weight over all files: ${totalMaxWeight}")
    } else {
        println("Maximum MC::Event weight over all files: none found")
    }
}