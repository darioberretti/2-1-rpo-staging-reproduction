function digest = sha256_file(filename)
%SHA256_FILE Return a lower-case SHA-256 digest without external programs.
%
% PROTOTYPE:
%   digest = rpo2ds.sha256_file(filename)
%
% DESCRIPTION:
%   Read the file in chunks using MATLAB/Java and close it on every exit.
%   Missing or unreadable inputs are recorded with text sentinels.
%
% INPUT:
%   filename  Input file path.
%
% OUTPUT:
%   digest    Lower-case hexadecimal string, or missing / unreadable.

if ~isfile(filename)
    digest = "missing";
    return
end
file_identifier = fopen(filename, 'r');
if file_identifier < 0
    digest = "unreadable";
    return
end
cleanup = onCleanup(@() fclose(file_identifier));
engine = java.security.MessageDigest.getInstance('SHA-256');
while ~feof(file_identifier)
    engine.update(fread(file_identifier, 8*1024*1024, '*uint8'));
end
digest = lower(string(reshape(dec2hex(typecast(engine.digest(), 'uint8'))', 1, [])));
end
