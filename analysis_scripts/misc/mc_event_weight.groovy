#!/usr/bin/env groovy

import org.jlab.jnp.hipo4.io.HipoReader
import org.jlab.jnp.hipo4.data.Event
import org.jlab.jnp.hipo4.data.Bank
import org.jlab.jnp.hipo4.data.SchemaFactory

if (args.length < 1) {
    System.err.println("Usage: groovy sum_mc_event_weights.groovy input.hipo")
    System.exit(1)
}

String inputFile = args[0]

HipoReader reader = new HipoReader()
reader.open(inputFile)

SchemaFactory factory = reader.getSchemaFactory()

if (!factory.hasSchema("MC::Event")) {
    System.err.println("ERROR: input file does not contain bank MC::Event")
    reader.close()
    System.exit(2)
}

Bank mcEvent = new Bank(factory.getSchema("MC::Event"))
Event event = new Event()

long nEvents = 0L
long nEventsWithMcEvent = 0L
long nRows = 0L
double sumWeights = 0.0

while (reader.hasNext()) {
    reader.nextEvent(event)
    event.read(mcEvent)

    nEvents++

    int rows = mcEvent.getRows()
    if (rows <= 0) {
        continue
    }

    nEventsWithMcEvent++

    for (int row = 0; row < rows; row++) {
        sumWeights += mcEvent.getFloat("weight", row)
        nRows++
    }
}

reader.close()

println("Input file: ${inputFile}")
println("Total events read: ${nEvents}")
println("Events with MC::Event: ${nEventsWithMcEvent}")
println("Total MC::Event rows summed: ${nRows}")
println("Sum of MC::Event weight: ${sumWeights}")