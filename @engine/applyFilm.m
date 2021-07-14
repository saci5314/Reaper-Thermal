function eng = applyFilm(eng, dm, injLoc)
    %%% This function models liquid or gaseous film cooling 

    eng.film = true;
    
    eng.dm_f = dm;
    eng.l_fInj = injLoc;
end