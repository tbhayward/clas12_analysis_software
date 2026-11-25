package exclusivepion;
/** Reject metadata/invalid run records without hiding real mixed-run input. */
public final class PhysicsRecord {
 public static boolean acceptRun(int run,int expectedRun){
  if(run<=0)return false;
  if(run!=expectedRun)throw new IllegalStateException("Input inventory run "+expectedRun+" does not match physics-event HIPO run "+run);
  return true;
 }
}
