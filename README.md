# TechnoSTM
## Un microscope atomique à effet tunnel

Le programme du Rasberry Pi 3b devra implémenter les fonctionnalités suivantes:

+ Contrôle du CAN ADS1115 pour l'acquisition du signal relatif à l'axe z

+ Contrôle du CNA AD5764 contrôler par un Arduino Uno pour le déplacement des actionneurs piézoélectriques : [Code Source](https://github.com/AlexandreBlais42/Arduino-EVAL-AD5764EB-pour-stm)

+ Acquisition et traitement de l'image 3D

+ 2 modes fonctionnement, 
    1. déplacement à hauteur constante
    2. déplacement à courant constant

+ Contrôle du moteur pas à pas pour la première approche grossière;

## Crédits

Un grand merci à [Dan Berard](https://dberard.com) pour son site qui a été d'une grande aide lors de la conception du microscope.
