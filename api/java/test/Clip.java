import rsf.RSF;
import rsf.Input;
import rsf.Output;

/* A simple Java program to clip a dataset.

Presently, there is no automatic self-documentation generation for use
with sfdoc.  Javadoc may be a better way to generate self-doc for Java
programs.

*/

public class Clip {
    static {
        System.loadLibrary("jrsf");
    }
    public static void main(String[] args){
        // Initialize command line argument passing
         RSF par = new RSF(args);
         // Get the input file name.
         Input input = new Input("in");
         Output output = new Output("out");
        // Get the value to clip to.
         float clip = par.getFloat("clip",0.0f);
        // Read our binary data.
        int n3 = input.getN(3);
        int n2 = input.getN(2);
        int n1 = input.getN(1);
        //Perform clipping operation on a single trace and write out.
        float[] data = new float[n1];
         for(int i = 0; i < n3; ++i){
            for(int j = 0; j < n2; ++j){
                input.read(data);
                for(int k = 0; k < n1; ++k){
                    if (data[k] > clip) data[k] = clip;
                    else if (data[k] < -clip) data[k] = -clip;
                }
                output.write(data);
            }
         }
         input.close();
         output.close();
    }
}
