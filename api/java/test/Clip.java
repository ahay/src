import rsf.Par;
import rsf.Reader;
import rsf.Writer;
import rsf.Header;

/* A simple Java program to clip a dataset.

Presently, there is no automatic self-documentation generation for use
with sfdoc.  Javadoc may be a better way to generate self-doc for Java
programs.

*/

public class Clip {
    public static void main(String[] args){
        // Initialize command line argument passing
         Par par = new Par(args);
         // Get the input file name.
         String input = par.getString("in","");
         // If the input file name is nothing, then quit!
         if (input.equals("")){
                System.out.println("Did not find input file!");
                System.exit(1);
         }
         //If the output file name is nothing, then quit!
         String output = par.getString("out","");
         if (output.equals("")){
                System.out.println("Did not find output file!");
                System.exit(1);
         }
        // Get the value to clip to.
         float clip = par.getFloat("clip",0.0f);
        //Read our header file.
         Header header = Reader.readHeader(input);
        // Read our binary data.
         float[][][] data = Reader.readBinary3D(header);
        //Initialize our array values.
         int n3 = header.n[2];
         int n2 = header.n[1];
         int n1 = header.n[0];
        //Perform clipping operation.
         for(int i = 0; i < n3; ++i){
            for(int j = 0; j < n2; ++j){
                for(int k = 0; k < n1; ++k){
                    float trace = data[i][j][k];
                    if (trace > clip) data[i][j][k] = clip;
                    else if (trace < -clip) data[i][j][k] = -clip;
                }
            }
         }
        //Write our data out, using the same header values to the file
        //located at: output.
         Writer.writeRSF(header,data,output);
    }
}
