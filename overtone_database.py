import os

flights = ['B737', 'B738', 'B739', 'B77W', 'B772', 'B788', 'B789', 'B763', 'B744','B733','B732','B77L','B748','CRJ2', 'A332', 'A359', 'A359', 'E75S','B190','BE20','C208','DH8A',
            'AT73','SW4','PC12','DH3T','C441','B18T','B350','BE10','AS50','R44', 'PA31','DHC2','GA8','C180','C182','C206','C172','PA32','PA46','CH7B','PA30','C46','BE35','PA18','PA34','C185']

for eq in  flights:
    file = 'output/Inversion_Results/'+eq+'data_atmosphere_full.txt'
   