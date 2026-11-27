package exclusivepion;
public final class PhysicsRecordTest {
 public static void main(String[] args){
  int accepted=0;
  for(int run:new int[]{0,0,6265,6265})if(PhysicsRecord.acceptRun(run,6265))accepted++;
  if(accepted!=2)throw new AssertionError("Metadata must not hide following data");
  if(PhysicsRecord.acceptRun(-1,6265))throw new AssertionError("Invalid run");
  try{PhysicsRecord.acceptRun(6266,6265);throw new AssertionError("Must reject real mixed-run input");}
  catch(IllegalStateException expected){}
  System.out.println("Record regression passed: repeated run-zero metadata, subsequent data, and real run mismatch");
 }
}
