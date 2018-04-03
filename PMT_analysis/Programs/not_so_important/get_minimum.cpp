void get_minimum(){
	fstream infile;
	infile.open( "Auswertung/a_o_Q1_T=+23C_007.dat" ) ;
	
	if ( infile.is_open() ) {
		cout << "File is open" << endl;
	}
	else{
		cout << "Failed to open file" << endl;
	}
	
	string line;
	float a = 1.0;
	float b = 1.0;
	float c = 1.0;
	float d = 1.0;
	float alpha = 1.0;
	float omega = 1.0;
	float Q1 = 1.0;
	float std_dev_Q1 = 1.0;
	while (getline(infile, line))
	{
		infile >> a >> b >> c >> d;
		if (d < std_dev_Q1){
			alpha = a;
			omega = b;
			Q1 = c;
			std_dev_Q1 = d;
		}
	}	
	cout << "<<<<< minimum found: <<<<<" << endl;
	cout << "alpha:\t\t" << alpha << endl;
	cout << "omega:\t\t" << omega << endl;
	cout << "std(Q1):\t" << std_dev_Q1 << endl;
	infile.close();	
}