import json
import sys
import urllib

# Credentials and feature service information
username = "koen.reef@hydrologic.com"
password = "yu7!JEc5M77#"
URL = "<the services url you want to connect to>"

# obtain a token
referer = "http://www.arcgis.com/"
query_dict = {"username": username, "password": password, "referer": referer}

query_string = urllib.parse.urlencode(query_dict)
url = "https://www.arcgis.com/sharing/rest/generateToken"
token = json.loads(urllib.urlopen(url + "?f=json", query_string).read())

if "token" not in token:
    print(token["error"])
    sys.exit(1)

query_dict = {"f": "json", "token": token["token"]}

# query your services url
jsonResponse = urllib.urlopen(URL, urllib.urlencode(query_dict))
