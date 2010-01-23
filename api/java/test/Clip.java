import rsf.Par;
import rsf.Reader;
import rsf.Writer;
import rsf.Header;

public class Clip {
    public static void main(String[] args){
         Par par = new Par(args);

         float clip = par.getFloat("clip",0.0f);
        
         Header header = Reader.readHeader("in");

         float[][][] data = Reader.readBinary3D(header);

         int n3 = header.n[2];
         int n2 = header.n[1];
         int n1 = header.n[0];

         for(int i = 0; i < n3; ++i){
            for(int j = 0; j < n2; ++j){
                for(int k = 0; k < n1; ++k){
                    float trace = data[i][j][k];
                    if (trace > clip) data[i][j][k] = clip;
                    else if (trace < -clip) data[i][j][k] = -clip;
                }
            }
         }

         Writer.writeRSF(header,data,"out");
    }
}
