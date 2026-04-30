# Déployer le site PT sur GitHub Pages avec `www.persistencetheory.org`

Ce guide explique comment publier le site Astro `PT-WEBSITE` sur GitHub Pages et le relier au domaine `www.persistencetheory.org`.

Objectif recommandé ici : utiliser `https://www.persistencetheory.org` comme domaine canonique, et laisser `https://persistencetheory.org` rediriger automatiquement vers `www`.

> État actuel du dépôt : le site est déjà configuré pour GitHub Pages en sortie statique, mais `astro.config.mjs`, `public/CNAME` et quelques textes internes pointent actuellement vers `https://persistencetheory.org` sans `www`. Si tu veux que `www.persistencetheory.org` soit le domaine officiel, applique les ajustements de l’étape 2.

## 1. Préparer le dépôt GitHub

Créer un dépôt GitHub, par exemple :

```text
persistencetheory.org
```

ou :

```text
PT-WEBSITE
```

Puis, depuis le dossier local du site :

```bash
cd "/Volumes/PT-YS-0326/LA THEORIE DE LA PERSITANCE/PT-WEBSITE"
git remote add origin git@github.com:TON_COMPTE_GITHUB/NOM_DU_DEPOT.git
git branch -M main
git push -u origin main
```

Si le remote existe déjà :

```bash
git remote -v
git remote set-url origin git@github.com:TON_COMPTE_GITHUB/NOM_DU_DEPOT.git
git push -u origin main
```

## 2. Configurer le domaine canonique dans le code

Pour publier officiellement sur `www.persistencetheory.org`, il faut que le site généré connaisse cette URL.

Dans `astro.config.mjs`, mettre :

```js
site: 'https://www.persistencetheory.org',
```

Dans `public/CNAME`, mettre exactement :

```text
www.persistencetheory.org
```

Vérifier les autres occurrences textuelles :

```bash
rg -n "persistencetheory.org"
```

Les occurrences importantes à harmoniser sont généralement :

```text
astro.config.mjs
public/CNAME
src/layouts/BaseLayout.astro
src/pages/about.astro
src/pages/en/about.astro
src/pages/monograph.astro
src/pages/en/monograph.astro
README.md
```

Après modification :

```bash
npm run build
git add astro.config.mjs public/CNAME src README.md
git commit -m "Configure canonical www domain"
git push
```

## 3. Activer le workflow GitHub Pages

Le dépôt contient déjà un workflow prêt, mais désactivé :

```text
.github/workflows/deploy.yml.disabled
```

Pour l’activer :

```bash
mv .github/workflows/deploy.yml.disabled .github/workflows/deploy.yml
git add .github/workflows/deploy.yml
git commit -m "Enable GitHub Pages deployment"
git push
```

Le workflow actuel fait :

1. installation Node 20 ;
2. `npm ci` ;
3. `npm run build` ;
4. upload de `dist/` ;
5. déploiement via GitHub Pages.

Il faut que `package-lock.json` soit bien commité, ce qui est déjà le cas dans ce projet.

## 4. Configurer GitHub Pages côté GitHub

Dans GitHub :

1. ouvrir le dépôt ;
2. aller dans `Settings` ;
3. ouvrir `Pages` ;
4. dans `Build and deployment`, choisir `Source: GitHub Actions` ;
5. attendre qu’un premier run se termine dans l’onglet `Actions`.

Quand le workflow réussit, GitHub Pages publie le contenu statique généré depuis `dist/`.

## 5. Vérifier le domaine dans GitHub

GitHub recommande de vérifier le domaine au niveau du compte ou de l’organisation pour éviter les prises de contrôle de domaine.

Dans GitHub :

1. aller dans les réglages du compte ou de l’organisation ;
2. ouvrir `Pages` ;
3. ajouter le domaine `persistencetheory.org` ;
4. GitHub donne un enregistrement `TXT` à ajouter chez le registrar/DNS ;
5. attendre la propagation ;
6. cliquer sur `Verify`.

Commande de vérification, en remplaçant `TON_COMPTE_GITHUB` :

```bash
dig _github-pages-challenge-TON_COMPTE_GITHUB.persistencetheory.org +nostats +nocomments +nocmd TXT
```

Garder ce TXT ensuite : il protège aussi les sous-domaines immédiats comme `www.persistencetheory.org`.

## 6. Configurer le domaine personnalisé dans le dépôt

Dans le dépôt GitHub :

1. `Settings` ;
2. `Pages` ;
3. dans `Custom domain`, entrer :

```text
www.persistencetheory.org
```

4. cliquer sur `Save`.

GitHub peut afficher une attente de vérification DNS. C’est normal tant que les enregistrements DNS ne sont pas propagés.

## 7. Configurer les DNS chez le registrar

Chez le fournisseur DNS de `persistencetheory.org`, créer les enregistrements suivants.

### Sous-domaine `www`

Pour que `www.persistencetheory.org` pointe vers GitHub Pages :

| Type | Nom | Valeur |
|---|---|---|
| `CNAME` | `www` | `TON_COMPTE_GITHUB.github.io` |

Ne pas mettre le nom du dépôt dans la cible CNAME. La cible doit être le domaine GitHub Pages par défaut du compte ou de l’organisation.

### Domaine racine `persistencetheory.org`

Pour que `persistencetheory.org` redirige correctement vers `www.persistencetheory.org`, ajouter les `A` records GitHub Pages :

| Type | Nom | Valeur |
|---|---|---|
| `A` | `@` | `185.199.108.153` |
| `A` | `@` | `185.199.109.153` |
| `A` | `@` | `185.199.110.153` |
| `A` | `@` | `185.199.111.153` |

Optionnel, pour IPv6 :

| Type | Nom | Valeur |
|---|---|---|
| `AAAA` | `@` | `2606:50c0:8000::153` |
| `AAAA` | `@` | `2606:50c0:8001::153` |
| `AAAA` | `@` | `2606:50c0:8002::153` |
| `AAAA` | `@` | `2606:50c0:8003::153` |

Supprimer les anciens enregistrements contradictoires pour `@` ou `www`, notamment :

- redirections web opaques ;
- anciens `A` records ;
- ancien `CNAME www` vers un autre service ;
- wildcard `*.persistencetheory.org`.

GitHub déconseille les wildcards DNS pour Pages, car ils augmentent le risque de prise de contrôle de sous-domaines.

## 8. Vérifier la propagation DNS

La propagation peut être immédiate ou prendre jusqu’à 24 heures.

Vérifier `www` :

```bash
dig www.persistencetheory.org +nostats +nocomments +nocmd
```

Tu dois voir une chaîne de type :

```text
www.persistencetheory.org.  CNAME  TON_COMPTE_GITHUB.github.io.
```

Vérifier l’apex :

```bash
dig persistencetheory.org +noall +answer -t A
```

Tu dois voir les IP GitHub Pages :

```text
185.199.108.153
185.199.109.153
185.199.110.153
185.199.111.153
```

## 9. Activer HTTPS

Quand GitHub a validé le domaine :

1. retourner dans `Settings > Pages` ;
2. cocher `Enforce HTTPS`.

L’option peut mettre du temps à devenir disponible après la configuration DNS.

## 10. Tester le site publié

Tester les pages principales :

```bash
curl -I https://www.persistencetheory.org/
curl -I https://www.persistencetheory.org/en/
curl -I https://www.persistencetheory.org/tableau-periodique/
curl -I https://www.persistencetheory.org/affinites-electroniques/
```

Tester la redirection apex vers `www` :

```bash
curl -I https://persistencetheory.org/
```

Selon la propagation et la configuration GitHub, on doit obtenir une réponse HTTPS valide et, idéalement, une redirection vers `https://www.persistencetheory.org/`.

## 11. Déploiement quotidien

Une fois tout configuré :

```bash
npm run build
git status
git add .
git commit -m "Update site"
git push
```

Le push sur `main` déclenche automatiquement GitHub Actions et republie le site.

## 12. Dépannage rapide

### La page GitHub Pages affiche une 404

Vérifier :

- `Settings > Pages > Source` doit être `GitHub Actions` ;
- le workflow `Deploy site to GitHub Pages` doit être vert ;
- le fichier `.github/workflows/deploy.yml` doit exister ;
- le build doit produire `dist/index.html`.

### Le domaine custom ne se valide pas

Vérifier :

- `public/CNAME` contient exactement `www.persistencetheory.org` ;
- le champ `Custom domain` dans GitHub Pages contient exactement `www.persistencetheory.org` ;
- le DNS `CNAME www` pointe vers `TON_COMPTE_GITHUB.github.io` ;
- aucun wildcard `*.persistencetheory.org` ne détourne la résolution.

### HTTPS n’est pas disponible

Attendre la propagation DNS, puis vérifier :

```bash
dig www.persistencetheory.org +nostats +nocomments +nocmd
dig persistencetheory.org +noall +answer -t A
```

Si les DNS sont bons, retourner dans `Settings > Pages` et réessayer `Enforce HTTPS`.

### Les liens ou le sitemap pointent vers l’apex sans `www`

Rechercher les URLs codées en dur :

```bash
rg -n "https://persistencetheory.org|persistencetheory.org"
```

Puis harmoniser vers :

```text
https://www.persistencetheory.org
```

Relancer :

```bash
npm run build
git add .
git commit -m "Normalize site canonical URL"
git push
```

## Sources officielles

- GitHub Docs, configuration de la source GitHub Pages : https://docs.github.com/en/pages/getting-started-with-github-pages/configuring-a-publishing-source-for-your-github-pages-site
- GitHub Docs, domaines personnalisés, DNS, `www`, apex et HTTPS : https://docs.github.com/en/pages/configuring-a-custom-domain-for-your-github-pages-site/managing-a-custom-domain-for-your-github-pages-site
- GitHub Docs, vérification du domaine pour éviter les prises de contrôle : https://docs.github.com/en/pages/configuring-a-custom-domain-for-your-github-pages-site/verifying-your-custom-domain-for-github-pages
- Astro Docs, déploiement Astro sur GitHub Pages : https://docs.astro.build/fr/guides/deploy/github/
