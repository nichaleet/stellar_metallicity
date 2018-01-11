function wc, file, lines=lines, words=words, characters=characters
    if ~keyword_set(lines) and ~keyword_set(words) and ~keyword_set(characters) then lines = 1
    spawn, 'wc '+file, wcout
    if keyword_set(lines) then return, long(strmid(wcout, 0, 6))
    if keyword_set(words) then return, long(strmid(wcout, 6, 7))
    if keyword_set(characters) then return, long(strmid(wcout, 13, 7))
end
