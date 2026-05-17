#include "common.h"


namespace orb {

    class OrbRecorder {

        // *******
        // @brief data recorder
        // *******

        public:
            OrbRecorder();
            ~OrbRecorder();

            void record_cartesian(timeposvel &posvel_data, const std::string &filename);
            void record_oelem(timeoelem &oelem_data, const std::string &filename);
            void record_geqoe(timegeqoe &geqoe_data, const std::string &filename);
    };




}