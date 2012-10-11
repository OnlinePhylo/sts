//
// Processing.org app for drawing a series of phyloSMC forest particles
//
int generation = 0;
int last_generation = 29;
int gen_time = 10;
int cur_time = gen_time+1;

float branch_alpha = 128;
float branch_shade = 255;

BufferedReader reader;
String filename = "/home/koadman/git/smctc/examples/phylo/viz_data.csv";
String frame_path = "/home/koadman/git/smctc/examples/phylo/movie/";
PFont f;

void setup() 
{
  size(640, 360);
  frameRate(30);
  try{
    reader = new BufferedReader (new FileReader (filename)) ; 
  }catch(Exception e){
    println("ERROR >> "+e);
  }
  f = createFont("Arial",16,true);
}

final int LEFT = 0;
final int RIGHT = 1;
final int LFROM = 2;
final int RFROM = 3;
final int TO = 4;

Vector< float[] > vizdata = new Vector< float[] >();

void draw() 
{
  // load data if needed
  cur_time++;
  if(cur_time > gen_time && generation < last_generation){
    cur_time = 0;
    generation++;
    vizdata = new Vector< float[] >();
  try {
    while (reader.ready()) {
      String line=reader.readLine();
      if(line.startsWith("#")){
         break;
      }
      StringTokenizer tok = new StringTokenizer(line);
      float[] data = new float[5];
      int i=0;
      while(tok.hasMoreTokens()){
        data[i] = Float.parseFloat(tok.nextToken());
        i++;
      }
      vizdata.addElement(data);
    }
  } catch (Exception e) {
    println("ERROR >> "+e);
  }

  }
  background(51);

  fill(255);
  textFont(f,12);
  text("Phylogenetic inference via Sequential Monte Carlo, 30 taxa, 1000 sites", 25, 30);
  textFont(f,10);
  text("Generation: ",25,55);
  text((new Integer(generation)).toString(),85,55);

  stroke(branch_shade, branch_alpha);
  for(int i = 0; i < vizdata.size(); i++){
    float[] data = vizdata.get(i);
      line(data[LEFT], height-data[LFROM], data[LEFT], height-data[TO]);
      line(data[RIGHT], height-data[RFROM], data[RIGHT], height-data[TO]);
      line(data[LEFT], height-data[TO], data[RIGHT], height-data[TO]);
  }
  String frame_name = frame_path + "frame_######.png";
  saveFrame(frame_name);
}

