#include "ProjData.h"

ProjData::ProjData() {}

ProjData::~ProjData() {}

bool ProjData::readSTIRProjData(const std::string& path) 
{
  header_path = path;
  ifstream input_file(path.c_str(), ios::in);
  
  string line;
  string sep = "=";
  string sep_comment = "#";
  string sep_elt = ","; 

  if (!input_file)
  {
    Cerr("***** castor-GATERootToCastor :: The input file is not read properly" << endl);
    Exit(EXIT_FAILURE);
  }
  while(getline(input_file, line))
  {
    if (line.find(sep_comment) != string::npos) 
      line = line.substr(0, line.find_first_of(sep_comment)) ; 

    //Sinogram binary file
    if (line.find("name of data file") != string::npos)
    {
      getSublineAfterTheSeparator(line, sep);
      binary_path = GetPathOfFile(header_path) + line;     
    }

    //Segments
    if (line.find("!matrix size [4]") != string::npos)
    {
      getSublineAfterTheSeparator(line, sep);
      if (ConvertFromString(line, &segments_no))
      {
        Cerr("***** castor-GATERootToCastor :: Exception when trying to value '" << line << " for option: segments number" << endl);
        Exit(EXIT_FAILURE);
      }

      int temp_segments_no = (segments_no+1)/2;
      for (int i=1; i<=temp_segments_no; i++)
        segments_sizes.push_back(i);
      for (int i=temp_segments_no-1; i>=1; i--)
        segments_sizes.push_back(i);
    }

    //Views
    if (line.find("!matrix size [3]") != string::npos)
    {
     getSublineAfterTheSeparator(line, sep); 
     if (ConvertFromString(line, &views_no))
      {
        Cerr("***** castor-GATERootToCastor :: Exception when trying to value '" << line << " for option: views number" << endl);
        Exit(EXIT_FAILURE);
      }
    }  

    //Tangential coordinates
    if (line.find("!matrix size [1]") != string::npos)
    {
      getSublineAfterTheSeparator(line, sep);
      if (ConvertFromString(line, &tangential_coord_no))
      {
        Cerr("***** castor-GATERootToCastor :: Exception when trying to value '" << line << " for option: tangential coordinates number" << endl);
        Exit(EXIT_FAILURE);
      } 
    }

    // Rings number
    if (line.find("Number of rings") != string::npos)
    {
      getSublineAfterTheSeparator(line, sep);
      if (ConvertFromString(line, &rings_no))
      {
        Cerr("***** castor-GATERootToCastor :: Exception when trying to value '" << line << " for option: rings number" << endl);
        Exit(EXIT_FAILURE);
      }
    }

    // Distance between rings
    if (line.find("Distance between rings (cm)") != string::npos)
    {
      getSublineAfterTheSeparator(line, sep);
      if (ConvertFromString(line, &distance_between_rings))
      {
        Cerr("***** castor-GATERootToCastor :: Exception when trying to value '" << line << " for option: distance between rings" << endl);
        Exit(EXIT_FAILURE);
      }
    }
      
    // Effective bin size
    if (line.find("effective central bin size (cm)") != string::npos)
    {
      getSublineAfterTheSeparator(line, sep);
      if (ConvertFromString(line, &effective_bin_size))
      {
        Cerr("***** castor-GATERootToCastor :: Exception when trying to value '" << line << " for option: effective bin size" << endl);
        Exit(EXIT_FAILURE);
      }
    }   
    
    axial_fov = distance_between_rings*static_cast<float>(rings_no);
    tangential_fov = effective_bin_size * static_cast<float>(tangential_coord_no);
  }
  // Print info about the data
  printInfoAboutProjData();

  // Read projection data
  readSTIRProjBinaryData();

  return true;
}



void ProjData:: printInfoAboutProjData()
{
  Cout(endl << "STIR SSS basic information abpout the projection data:"<< endl);
  Cout("Segments number: " << segments_no << endl);
  Cout("Segments sizes: {");Cout("Segments sizes: {");
  for (int i=0; i<segments_sizes.size(); i++)
    Cout(segments_sizes[i] << " ");
  Cout("}" << endl);
  Cout("Views number: " << views_no << endl);
  Cout("Tangential coordinates number: " << tangential_coord_no << endl);
  Cout("Tangential FOV: "<<tangential_fov<<" cm"<<endl);
  Cout("Rings number: " << rings_no << endl);
  Cout("Distance between rings: " << distance_between_rings << " cm"<< endl);
  Cout("Axial FOV: " << axial_fov << " cm" << endl);
  Cout("Effective bin size: " << effective_bin_size << " cm" << endl);
  Cout("Binary sinogram file: " << binary_path << endl);

}



void ProjData::readSTIRProjBinaryData()
{
  ifstream proj_file(binary_path.c_str(), std::ios::binary | std::ios::in);
  std::vector<float> binary_data;
  
  float f;

  if (proj_file)
  {
    while (proj_file.read((char*)&f, sizeof(float)))
    {
      binary_data.push_back(f);
    }
    Cout("Number of bins: " << binary_data.size() << endl); 
  }
  else
  {
    Cerr("***** castor-GATERootToCastor :: Exception when trying to open STIR projection binary file: " << binary_path << endl);
    Exit(EXIT_FAILURE); 
  }

  //Check the data consistency
  int sum_of_segments_elements = std::accumulate(segments_sizes.begin(), segments_sizes.end(), 0);
  Cout ("Sum of segments elements: " << sum_of_segments_elements << endl);
  if(!tangential_coord_no*sum_of_segments_elements*views_no == binary_data.size())
  {
    Cerr("***** castor-GATERootToCastor :: The projection data (binary file size) is inconsistent with the header information!" << endl);
    Exit(EXIT_FAILURE);
  }
   
  // Save data to segments
  uint32_t bins_counter = 0;
  for (int s=0; s<segments_sizes.size(); s++)
  {
    std::vector<std::vector<std::vector<float>>> segment_data(tangential_coord_no, std::vector<std::vector<float>>(segments_sizes[s], std::vector<float>(views_no, 0.0)));
    
    for (int k=0; k<views_no; k++)
    {
      for (int j=0; j<segments_sizes[s]; j++)
      {
        for(int i=0; i<tangential_coord_no; i++)
        { 
          segment_data[i][j][k] = binary_data[bins_counter]; 
          bins_counter++;
        }
      }
    }
    //Save data to segments properly
    segments.push_back(Segment(segment_data));
  }
}



float ProjData::getBinValue(uint32_t segment_id, uint32_t sinogram_id, uint32_t view_id, uint32_t tangential_id)
{
  return segments[segment_id].getBinValue(tangential_id, sinogram_id, view_id);
}



void ProjData::getSublineAfterTheSeparator(std::string &line, std::string &sep)
{
  line = line.substr(line.find_first_of(sep)+1);
  line.erase(0, line.find_first_not_of(" !\t\r\n")); // Erase all blank stuff before the first character
}



float ProjData::getScatterNumber(float posX1, float posY1, float posZ1, float posX2, float posY2, float posZ2)
{
  uint32_t segment_id, sinogram_id, view_id, tangential_id, ringID1, ringID2;
  float angle, line_to_center_point_distance, directional_parameter;
  // Overcome the issue with the coincidences at the boundaries of the scanner FOV (wrong ring ID is assigned)
  if (posZ1/10. >= axial_fov/2)
    posZ1 = axial_fov*10./2 - distance_between_rings*10./100.;   // convert half of the axial fov to mm (*10) and subtract 1% of the distance between rings 
  if (posZ2/10. >= axial_fov/2)
    posZ2 = axial_fov*10./2 - distance_between_rings*10./100.;   // convert half of the axial fov to mm (*10) and subtract 1% of the distance between rings 
  if (posZ1/10. <= -axial_fov/2)
    posZ1 = -axial_fov*10./2 + distance_between_rings*10./100.;   // convert half of the axial fov to mm (*10) and add 1% of the distance between rings 
  if (posZ2/10. <= -axial_fov/2)
    posZ2 = -axial_fov*10./2 + distance_between_rings*10./100.;   // convert half of the axial fov to mm (*10) and add 1% of the distance between rings 

  // Calculate the ringsIDs based on the axial fov and distance between rings information
  ringID1 = static_cast<int>(posZ1/10. + axial_fov/2) / distance_between_rings;
  ringID2 = static_cast<int>(posZ2/10. + axial_fov/2) / distance_between_rings;

  // segment_id
  segment_id = (ringID1 - ringID2) + rings_no - 1;
  
  // sinogram_id
  // The sinogram_id will be smaller value from the ringIDs
  sinogram_id = std::min(ringID1, ringID2);   

  calculatePointLineDistance(posX1/10., posY1/10., posX2/10., posY2/10., line_to_center_point_distance, directional_parameter);
  angle = get2DAngle(posX1/10., posY1/10., posX2/10., posY2/10.);
  if (directional_parameter<0.)
    angle = 180.0 -angle;

  // view_id
  view_id = static_cast<int>(angle / (180./views_no));

  // tangential_id
  if (line_to_center_point_distance < -tangential_fov/2.)
    tangential_id = tangential_coord_no-1;
  else if (line_to_center_point_distance >= tangential_fov/2.)
    tangential_id = 0;
  else
    tangential_id = static_cast<int> ((tangential_fov/2.-line_to_center_point_distance) / (effective_bin_size));
  
  return getBinValue(segment_id, sinogram_id, view_id, tangential_id); 
}



float ProjData::get2DAngle(float posX1, float posY1, float posX2, float posY2)
{
  // emission point is a point to build the triangle for the cosine law with the coordinates (posX1, posY2)
  float dist_ep_c1 = abs(posY2-posY1); // distance emission - registration 1
  float dist_ep_c2 = abs(posX2-posX1); // distance emission - registration 2
  float dist_c1_c2 = pow((posX2-posX1)*(posX2-posX1) + (posY2-posY1)*(posY2-posY1), 0.5); // distance registration 1 - registration 2
  // Cosine law
  float cosBeta = (dist_ep_c1*dist_ep_c1 + dist_c1_c2*dist_c1_c2 - dist_ep_c2*dist_ep_c2)/(2*dist_ep_c1*dist_c1_c2);

  if (cosBeta>=1.0)
    cosBeta = 0.9999;  //not 1 or -1 to overcome the issue once the 180 degree angle is calculated and index of view_id will be greater than the size of the matrix 
  if (cosBeta<=-1.0)
    cosBeta = -0.9999; 
  return acos(cosBeta)*180.0/M_PI;
}



void ProjData::calculatePointLineDistance(float x1, float y1, float x2, float y2, float &line_to_center_point_distance, float &directional_parameter)
{

  //using the equation: d = |A*x_0 + B*y_0 + C| / sqrt(A^2+B^2)
  //where (x_0, y_0) is a beginning of the coordinate system (0, 0)
  //and the line between coincidence points is given as y=ax+b and 
  //transformed to the A*x+B*y+C=0. Assuming both:
  // d = |b| / sqrt(a*a+1)
  //where b and a are calculated from the system of two equations
  if (x1!=x2)
  {
    directional_parameter = (y1-y2)/(x1-x2);
    float b = y1-x1*(y1-y2)/(x1-x2);
    // absolute value of b is not used as it will be an indicator for the bin position calculation
    line_to_center_point_distance = b/pow(directional_parameter*directional_parameter+1, 0.5);
  }
  else
  {
   directional_parameter = 1;
   line_to_center_point_distance = x1;
  }
}
