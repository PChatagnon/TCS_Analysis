#ifndef Table_Class
#define Table_Class

using namespace std;

class Latex_Table_writter
{
public:
    string cs_latex_string = "";
    string e_cs_latex_string = "";
    string variable_latex_string = "";
    string e_variable_latex_string = "";

    string variable = "";

    string output_name = "";
    string output_folder = "";
    string output_content = "";

    Latex_Table_writter(string in_output_name, string in_output_folder, string in_variable)
    {
        output_name = in_output_name;
        output_folder = in_output_folder;
        variable = in_variable;
    }

    void Set_output_name(string in_output_name)
    {
        output_name = in_output_name;
    }

    void Set_output_folder(string in_output_folder)
    {
        output_folder = in_output_folder;
    }

    void Set_variable(string in_variable)
    {
        variable = in_variable;
    }

    void Add_value(int in_i, int in_size_array, double in_cs_value, double in_e_cs_value, double in_variable_value, double in_e_variable)
    {
        cs_latex_string += std::to_string(in_cs_value);
        e_cs_latex_string += std::to_string(in_e_cs_value);
        variable_latex_string += std::to_string(in_variable_value);
        e_variable_latex_string += std::to_string(in_e_variable);
        if (in_i != (in_size_array - 1))
        {
            cs_latex_string += " & ";
            e_cs_latex_string += " & ";
            variable_latex_string += " & ";
            e_variable_latex_string += " & ";
        }
    }

    void Format()
    {
        output_content = variable_latex_string + "\n" +
                         e_variable_latex_string + "\n" +
                         cs_latex_string + "\n" +
                         e_cs_latex_string + "\n";
    }

    void Save()
    {
        std::ofstream outputFile(output_folder + "/" + output_name);
        if (outputFile.is_open())
        {
            outputFile << output_content;
            outputFile.close();
            std::cout << "Table of the results saved to file: " << output_name << std::endl;
        }
        else
        {
            std::cerr << "Error opening file for writing table at "<< output_name << std::endl;
        }
    }
};

#endif