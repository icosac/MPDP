#include "DPS.hh"

int plot_dps(const std::string& filename){
		std::ifstream file(filename);
		if (!file.is_open()) {
			std::cout << "Error opening file" << std::endl;
			return 1;
		}

		int n_lines = 0;
		file >> n_lines;

		std::vector<Configuration2> points;

		double x, y;

		while(file >> x >> y){
			points.push_back (Configuration2 (x, y, ANGLE::FREE));
		}

		for (size_t i = 0; i < points.size() - 1; i++){
			points[i].th(atan2( points[i+1].y() - points[i].y(), points[i+1].x() - points[i].x()));
		}
		points.back().th(points[points.size()-2].th());

		K_T kmax = 1.0;

		std::cout << std::endl;
		std::ofstream draw_file("Dubins3PSquares.asy");
		for (size_t i = 0; i < points.size()-1; i++){
			Dubins d(points[i], points[i+1], {kmax});
			std::cout << "=============\n" << d << std::endl;
			d.draw(draw_file, std::to_string(i), 500, 500, false, false, i==0);
		}
		draw_file.close();

		return 0;
}