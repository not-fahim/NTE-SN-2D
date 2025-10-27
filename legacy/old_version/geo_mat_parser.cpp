#include "geo_mat_parser.hpp"

using namespace std;

void parseCoordinates( stringstream& ss,  vector<cell_class>& targetVector) {
    double x;
    int y;
    char openParen, delimiter, closeParen;
    cell_class cell;

    while (ss >> openParen && openParen == '(') 
    {
        ss >> x; // Read the first number (x-coordinate)

        //see if it's a comma
        if (ss.peek() == ',') {
            ss >> delimiter; // read and move ss to next
        }

        ss >> y; //read the meshnum and move to next
        ss >> closeParen; //read the closing paranteses and move to the next

        if (closeParen == ')') 
        {
            cell.cell_dim = x;
            cell.mesh_num = y;
            targetVector.push_back(cell);
        }
    }
}

std::string read_global_file()
{
    string filename = "input.dotb";
    ifstream inputfile(filename);
    string line;
    if(!inputfile.is_open())
    {
        cerr<<"error: Could not open global input file"<<endl;
        return "nan";
    }

    while(getline(inputfile, line))
    {
        if(line.empty())
            continue;
        
        stringstream linestream(line);
        string keyword;
        linestream >> keyword;

        if(keyword[0] == '%')
            continue;
        
        else if(keyword == "problem")
        {
            string problemfilename;
            linestream >> problemfilename;
            inputfile.close();
            return problemfilename;
        }
    }
    inputfile.close();
    return "nan";
}


input_class read_input_file()
{
    input_class input_obj;

    string problemfilename = read_global_file();
    ifstream inputFile(problemfilename); 
    
    if (!inputFile.is_open()) {
         cerr << "Error: Could not open:"<<problemfilename <<  endl;
        input_obj.name = "nan";
        return input_obj;
    }

    string line;
    // one line at a time
    while ( getline(inputFile, line)) 
    {
        // Skip any empty lines
        if (line.empty()) 
        {
            continue; 
        }
        stringstream linestream(line);
        string keyword;
        linestream >> keyword; //first word to identify the line's purpose

        if(keyword[0] == '%') //comment line
        {
            continue;
        }
        // the line based on its keyword
        if (keyword == "Name") 
        {
            linestream >> input_obj.name;
        }
        else if(keyword=="tol")
        {
            double value;
            if(linestream >> value)
                input_obj.tol_in = value;
            if(linestream >> value)
                input_obj.tol_out = value;
        }
        else if(keyword=="maxit")
        {
            double value;
            if(linestream >> value)
                input_obj.max_it = value;
        }
        else if (keyword == "CellX") 
        {
            parseCoordinates(linestream, input_obj.cellX_vector);
        } 
        else if (keyword == "CellY") 
        {
            parseCoordinates(linestream, input_obj.cellY_vector);
        } 
        else if(keyword == "refinement")
        {
            int value;
            if(linestream >> value)
                input_obj.refinement = value;
        }
        else if(keyword == "BC")
        {
            double value;
            while(linestream >> value)
                input_obj.BC.push_back(value);
        }
        else if(keyword == "SN")
        {
            linestream >> input_obj.S_n;
        }
        else if(keyword == "xsfile")
        {
            linestream >> input_obj.xsfile;
        }
        else if (keyword == "MatMap") 
        {
            // Once "MatMap" is found, read all subsequent lines as matrix rows
            // until we hit an empty line or the end of the file.
            string matrixLine;
            while ( getline(inputFile, matrixLine) && !matrixLine.empty()) 
            {
                stringstream matrixStream(matrixLine);
                vector<int> row;
                int value;
                while (matrixStream >> value) 
                {
                    row.push_back(value);
                }
                if (!row.empty()) 
                {
                    input_obj.matmap_cell.push_back(row);
                }
            }
            std::reverse(input_obj.matmap_cell.begin(), input_obj.matmap_cell.end());
        }
    }
    inputFile.close();
    input_obj.mat_vector = read_material_input(input_obj.xsfile);
    if(input_obj.mat_vector.size()==0)
    {
        cerr<<"failed to read material object";
        input_obj.name = "no_mat";
    }
    return input_obj;

}

std::vector<mat_class> read_material_input(std::string filename)
{
    vector<mat_class> materials;
    double value;
    ifstream inputfile(filename);

    if(!inputfile.is_open())
    {
        cerr<<"error: Could not open material input file"<<endl;
        return materials;
    }
    
    string line;
    while(getline(inputfile, line))
    {
        if(line.empty())
            continue;
        
        stringstream linestream(line);
        string keyword;
        linestream >> keyword;

        if(keyword[0] == '%')
            continue;
        
        else if(keyword == "mat")
        {
            mat_class mat;
            linestream >>keyword;
            mat.name=keyword;
            //found a material and now evaluate the XS
            string matline;
            while(getline(inputfile, matline) )
            {

                stringstream matstream(matline);
                matstream >> keyword;
                if(keyword[0] == '%' || matline.empty())
                    continue;

                else if(keyword == "tot")
                {
                    while(matstream >> value)
                    {
                        mat.sigma_t.push_back(value);
                    }
                }
                else if(keyword == "nuf")
                {
                    while(matstream >> value)
                    {
                        mat.nu_sigma_f.push_back(value);
                    }
                }
                else if(keyword == "chi")
                {
                    while(matstream >> value)
                    {
                        mat.chi.push_back(value);
                    }
                }
                else if(keyword == "sca")
                {
                    string matrixLine;
                    while (getline(inputfile, matrixLine) && !matrixLine.empty())
                    {
                        stringstream matrixStream(matrixLine);
                        vector<double> row;
                        while (matrixStream >> value) 
                        {
                            row.push_back(value);
                        }
                        if (!row.empty()) 
                        {
                            mat.sigma_s.push_back(row);
                        }
                    }
                    break;
                }
                
            }
            materials.push_back(mat);
        }
    }
    inputfile.close();
    return materials;
}

void print_input_data(const input_class& data) {
    if (data.name == "nan") {
        std::cout << "Failed to read input file." << std::endl;
        return;
    }
    if (data.name == "no_mat") {
        std::cout << "after reading geometry, Failed to read material input file." << std::endl;
        return;
    }

    std::cout << "--- Successfully Parsed Geometry Data ---" << std::endl;
    std::cout << "Name: " << data.name << std::endl;
    std::cout << "Material File: " << data.xsfile << std::endl;
    std::cout << "SN Value: " << data.S_n << std::endl;

    std::cout << "\nCellX Data:" << std::endl;
    for (const auto& cell : data.cellX_vector) {
        std::cout << "  Dimension: " << cell.cell_dim << ", Mesh Number: " << cell.mesh_num << std::endl;
    }

    std::cout << "\nCellY Data:" << std::endl;
    for (const auto& cell : data.cellY_vector) {
        std::cout << "  Dimension: " << cell.cell_dim << ", Mesh Number: " << cell.mesh_num << std::endl;
    }

    std::cout << "\nBoundary Conditions (BC): ";
    for (double val : data.BC) {
        std::cout << val << " ";
    }
    std::cout << std::endl;

    std::cout << "\nMaterial Map:" << std::endl;
    for (const auto& row : data.matmap_cell) {
        std::cout << "  ";
        for (int val : row) {
            std::cout << val << " ";
        }
        std::cout << std::endl;
    }
    std::cout << "---------------------------------------" << std::endl;

    // --- New section to print material data ---
    std::cout << "\n--- Successfully Parsed Material Data ---" << std::endl;
    std::cout << std::fixed << std::setprecision(8); // For better float output

    for (const auto& mat : data.mat_vector) {
        std::cout << "\n## Material: " << mat.name << " ##" << std::endl;

        std::cout << "Total (sigma_t):      ";
        for (double val : mat.sigma_t) std::cout << val << "  ";
        std::cout << std::endl;

        std::cout << "Nu-Fission (nu_sigma_f):";
        for (double val : mat.nu_sigma_f) std::cout << val << "  ";
        std::cout << std::endl;

        std::cout << "Fission Spectrum (chi): ";
        for (double val : mat.chi) std::cout << val << "  ";
        std::cout << std::endl;

        std::cout << "Scattering Matrix (sigma_s):" << std::endl;
        for (const auto& row : mat.sigma_s) {
            std::cout << "  ";
            for (double val : row) {
                std::cout << val << "  ";
            }
            std::cout << std::endl;
        }
    }
    std::cout << "---------------------------------------" << std::endl;
}

