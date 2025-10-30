#!/bin/bash

# --- Configuration ---
OVPN_PATH="/home/cemal/Downloads/CemalKoba.ovpn"
VPN_PASSWORD=
SSH_KEY="/home/cemal/cemal_key"
SSH_USER="CemalKoba"
SSH_HOST=

# --- Connect to VPN ---
echo "[+] Importing VPN config..."
openvpn3 config-import --config "$OVPN_PATH" --persistent > /dev/null 2>&1

echo "[+] Starting VPN session..."
# Use expect to automatically provide password
expect <<EOF
spawn openvpn3 session-start --config "$OVPN_PATH"
expect "Password:"
send "$VPN_PASSWORD\r"
expect eof
EOF

# --- Check VPN status ---
echo "[+] Checking VPN sessions..."
openvpn3 sessions-list

# --- Connect to remote server ---
echo "[+] Connecting via SSH..."
ssh -i "$SSH_KEY" "$SSH_USER@$SSH_HOST"
