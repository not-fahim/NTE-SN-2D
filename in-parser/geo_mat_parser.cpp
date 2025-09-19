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

input_class read_input_file()
{
    input_class input_obj;
    string filename;
    cout<<" give name for the input file"<<endl;
    cin>> filename;
    //openfile
    ifstream inputFile(filename); 
    if (!inputFile.is_open()) {
         cerr << "Error: Could not open input.txt" <<  endl;
        input_obj.name = "nan";
        return input_obj;
    }

    string line;
    // 2. Read the file one line at a time
    while ( getline(inputFile, line)) 
    {
        // Skip any empty lines
        if (line.empty()) 
        {
            continue; 
        }

        stringstream lineStream(line);
        string keyword;
        lineStream >> keyword; // Read the first word to identify the line's purpose
        if(keyword[0] == '%') //comment line
        {
            continue;
        }
        // 3. Process the line based on its keyword
        if (keyword == "Name") 
        {
            lineStream >> input_obj.name;
        } 
        else if (keyword == "CellX") 
        {
            parseCoordinates(lineStream, input_obj.cellX);
        } 
        else if (keyword == "CellY") 
        {
            parseCoordinates(lineStream, input_obj.cellY);
        } 
        else if(keyword == "BC")
        {
            double value;
            while(lineStream >> value)
                input_obj.BC.push_back(value);
        }
        else if(keyword == "SN")
        {
            lineStream >> input_obj.S_n;
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
        }
    }

    inputFile.close();
    return input_obj;

}