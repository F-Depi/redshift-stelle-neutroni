## Introduzione
 Progetto finale del corso di Fisica Computazionale dell'Università di Trento,
 aa 2024

## Utilizzo
Partire da ```script/plot.py```, lì per ogni immagine che si può generare è
indicato quale script in C genera i dati.

Gli script in C che generano i dati sono:
 - ```test.c``` per verificare che il codice funzioni e le equazioni siano
 scritte correttamente. Controlla che il metodo numerico per invertire $P(\rho)$
 funzioni correttamente, confronta le energie delle 3 politropiche, verifica che
 gli step fatti da $\texttt{RK4}$ risolvano le TOV fino ad arrivare a $P \sim 0$
 e verifica la dipendenza del risultato dallo step di integrazione $h$.
 - ```main.c``` genera i dati per il grafico MR ($\texttt{get\\_MR()}$ è
 abbastanza lenta), genera tutti i dati $(r,~P(r),~m(r))$ per le stelle di massa
 massima e poi calcola il potenziale gravitazionale con quest'ultimi.
 - ```radianza.c``` calcola la radianza delle 3 stelle, controlla i giusti
 parametri da utilizzare per fare l'integrale relativo alla potenza e poi li usa
 per calcolare la potenza in funzione del raggio e la temperatura efficacie in
 funzione di quella propria della stella.
 - ```punto7.c``` calcola il grafico della temperatura efficacie rispetto alla
 pressione centrale della stella, usa i dati generati da ```main.c``` e i
 parametri per fare l'integrale trovati in ```radianza.c```.

Il file con le funzioni che servono a risolvere le TOV è ```fun.c```, il
corrispettivo header ```fun.h``` contiene anche i valori delle costanti $R0$,
$M0$, $P0$, etc, utilizzate per rendere le variabili adimensionali.
