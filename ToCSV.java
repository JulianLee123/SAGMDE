import java.io.BufferedReader;
import java.io.FileNotFoundException;
import java.io.FileReader;
import java.io.IOException;
import java.io.File;
import java.io.PrintWriter;


public class ToCSV {
	public static void main(String[] args){
		//Args[0]: file name excluding .txt extension
		//Args[1]: number of iterations per algorithm tested
		String name = args[0];
		StringBuilder[] sb=new StringBuilder[150];
		for(int i = 0; i < sb.length; i++){
			sb[i] = new StringBuilder();
		}
		//Setup info
		String[] metrics = new String[] {"Time (ms)","NMI"};
		for(int i = 0; i < 2; i++){
			sb[i*(Integer.parseInt(args[1])+1)+1].append(metrics[i]);
		}
		File file = new File(name + ".txt");
		//process input file
		try {
			BufferedReader br = new BufferedReader(new FileReader(file));
			String st = br.readLine(), alg, dataset = "";
			while (st != null){
				if(st.length() > 0 && st.charAt(0) == '-'){
					//new column entry
					st = br.readLine();//read space
					st = br.readLine();//read the header
					int spaceLoc = st.indexOf(' ');
					alg = st.substring(0, spaceLoc);
					if(!st.substring(spaceLoc + 1,st.length()).equals(dataset)){
						//add new dataset header
						dataset = st.substring(spaceLoc + 1,st.length());
						sb[0].append(dataset);
						for(int i = 0; i < sb.length; i++){
							sb[i].append(",");
						}
					}
					while(!(st = br.readLine()).equals("Results:"));
					//read in the actual data
					st = br.readLine();
					sb[0].append(alg + ",");
					int i = 1;
					for(; st != null && (st.length() == 0 || st.charAt(0) != '-'); i++){
						if(st.length() == 0 || st.charAt(0) != '*'){
							sb[i].append(st + ",");
						}
						st = br.readLine();
					}
				}
				else{
					st = br.readLine();
				}
			}
			br.close();
		} catch (FileNotFoundException e) {
			e.printStackTrace();
		} catch (IOException e) {
			e.printStackTrace();
		}
		//write results to csv
		try (PrintWriter writer = new PrintWriter(new File(name + ".csv"))) {
			for(int i = 0; i < sb.length; i++){
				writer.write(sb[i].toString());
				writer.write("\n");
			}
		} catch (FileNotFoundException e) {
			System.out.println(e.getMessage());
		}
	}
}
